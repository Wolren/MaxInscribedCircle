//! Score the native solvers against GEOS across a corpus, with an independent
//! check on both answers.
//!
//! ```text
//! GEOS_LIB_DIR=<lib dir> GEOS_INCLUDE_DIR=<include dir> GEOS_VERSION=3.14.1 \
//!   cargo run --release --features bench-geos --example parity_geos -- [dir ...]
//! ```
//!
//! Inputs: every `.wkt` file in the given directories (defaults to the local
//! georust-geo and GEOS test fixtures), plus closed-form shapes and a
//! deterministic synthetic corpus.
//!
//! Both solvers get the same tolerance, [`grid_tolerance`], so the comparison
//! measures the solvers rather than the tolerance gap between them. Two
//! measurements are taken without either solver's help and applied to both:
//!
//! - clearance: the true distance from the returned center to the nearest
//!   boundary segment. A radius above it means the reported circle does not fit
//!   inside the polygon.
//! - floor: the largest clearance found on a dense grid of interior points.
//!   Every one of those points is a feasible circle center, so a floor above a
//!   solver's radius proves that solver left radius on the table.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;
use std::str::FromStr;
use std::time::Instant;

use geo::Contains;
use geo_types::{Coord, LineString, Point, Polygon};
use wkt::Wkt;

use max_inscribed_circle::mic::grid_tolerance;
use max_inscribed_circle::mic::index::NearestBoundaryIndex;
use max_inscribed_circle::mic::input::{HostPolygon, SegmentIndex};
use max_inscribed_circle::mic::solver::bench_geos::solve_mic_geos_with_tolerance;
use max_inscribed_circle::mic::{maximum_inscribed_circle, MicOptions, MicResult};

// Fixture directories are taken from the environment so the harness does not
// carry one machine's layout. Pass directories as arguments, or set
// MIC_WKT_DIRS to a `;`-separated (Windows) or `:`-separated list.
fn default_wkt_dirs() -> Vec<PathBuf> {
    match std::env::var("MIC_WKT_DIRS") {
        Ok(raw) if !raw.trim().is_empty() => std::env::split_paths(&raw).collect(),
        _ => Vec::new(),
    }
}

/// Components above this vertex count skip the sampling floor (it sweeps the
/// ring per probe and is not worth the wait on giants).
const FLOOR_VERTEX_LIMIT: usize = 5_000;
/// Components larger than this are skipped entirely.
const COMPONENT_VERTEX_LIMIT: usize = 60_000;
/// Cap on real-world components pulled from the WKT corpus.
const REAL_COMPONENT_LIMIT: usize = 1500;
const FLOOR_GRID: usize = 64;

struct Rng(u64);

impl Rng {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

struct Case {
    name: String,
    kind: &'static str,
    polygon: Polygon<f64>,
    exact_radius: Option<f64>,
}

struct Outcome {
    name: String,
    kind: &'static str,
    verts: usize,
    native: Option<MicResult>,
    geos: Option<MicResult>,
    native_err: Option<String>,
    geos_err: Option<String>,
    clearance_native: f64,
    clearance_geos: f64,
    inside_native: bool,
    inside_geos: bool,
    floor: Option<f64>,
    floor_proves_native: Option<f64>,
    floor_proves_geos: Option<f64>,
    engine_native: String,
    approach_native: String,
    engine_geos: String,
    candidates_native: usize,
    candidates_geos: usize,
    unclosed_rings: usize,
    us_native: f64,
    us_geos: f64,
    exact_radius: Option<f64>,
}

fn ring(points: &[(f64, f64)]) -> LineString<f64> {
    LineString::from(
        points
            .iter()
            .map(|&(x, y)| Coord { x, y })
            .collect::<Vec<_>>(),
    )
}

fn regular_polygon(n: usize, radius: f64) -> Polygon<f64> {
    let mut pts = Vec::with_capacity(n + 1);
    for i in 0..n {
        let a = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        pts.push((radius * a.cos(), radius * a.sin()));
    }
    pts.push(pts[0]);
    Polygon::new(ring(&pts), vec![])
}

fn analytic_cases() -> Vec<Case> {
    let mut out = Vec::new();
    let mut push = |name: &str, polygon: Polygon<f64>, r: f64| {
        out.push(Case {
            name: name.to_string(),
            kind: "analytic",
            polygon,
            exact_radius: Some(r),
        });
    };

    push(
        "square 10",
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        5.0,
    );
    push(
        "rect 8x2",
        Polygon::new(
            ring(&[(0.0, 0.0), (8.0, 0.0), (8.0, 2.0), (0.0, 2.0), (0.0, 0.0)]),
            vec![],
        ),
        1.0,
    );
    push(
        "equilateral tri side 10",
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (5.0, 8.660_254_037_844_387),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        10.0 / (2.0 * 3.0f64.sqrt()),
    );
    for n in [5usize, 8, 32, 128] {
        push(
            &format!("regular {n}-gon R10"),
            regular_polygon(n, 10.0),
            10.0 * (std::f64::consts::PI / n as f64).cos(),
        );
    }
    // L shape, arm width 2: tangent to the two outer edges and the reflex
    // vertex, r = w*sqrt(2) / (1 + sqrt(2)).
    push(
        "L shape arm 2",
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 2.0),
                (2.0, 2.0),
                (2.0, 6.0),
                (0.0, 6.0),
                (0.0, 0.0),
            ]),
            vec![],
        ),
        2.0 * 2.0f64.sqrt() / (1.0 + 2.0f64.sqrt()),
    );
    // Square frame 10 with a concentric 4x4 hole: the corner candidate wins,
    // r = h*sqrt(2) / (1 + sqrt(2)) with h = (S - s)/2 = 3.
    push(
        "frame 10 hole 4",
        Polygon::new(
            ring(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 10.0),
                (0.0, 10.0),
                (0.0, 0.0),
            ]),
            vec![ring(&[
                (3.0, 3.0),
                (7.0, 3.0),
                (7.0, 7.0),
                (3.0, 7.0),
                (3.0, 3.0),
            ])],
        ),
        3.0 * 2.0f64.sqrt() / (1.0 + 2.0f64.sqrt()),
    );

    out
}

fn convex_hull(points: &[(f64, f64)]) -> Vec<(f64, f64)> {
    fn cross(o: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    }
    if points.len() < 3 {
        return Vec::new();
    }
    let mut pts = points.to_vec();
    pts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut lower: Vec<(f64, f64)> = Vec::new();
    for &p in &pts {
        while lower.len() >= 2 && cross(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }
    let mut upper: Vec<(f64, f64)> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower.push(lower[0]);
    lower
}

fn synthetic_cases(count: usize) -> Vec<Case> {
    let mut out = Vec::new();
    let mut rng = Rng(0x5eed_1234_9abc_def0);

    // Star-shaped: sorted angles around a center, random radii. Simple by
    // construction, usually non-convex.
    for i in 0..count / 3 {
        let k = 5 + (i % 10);
        let mut angles: Vec<f64> = (0..k)
            .map(|_| rng.range(0.0, 2.0 * std::f64::consts::PI))
            .collect();
        angles.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut pts: Vec<(f64, f64)> = angles
            .iter()
            .map(|&a| {
                let r = rng.range(1.0, 4.0);
                (r * a.cos(), r * a.sin())
            })
            .collect();
        pts.push(pts[0]);
        out.push(Case {
            name: format!("star {i} k={k}"),
            kind: "synthetic",
            polygon: Polygon::new(ring(&pts), vec![]),
            exact_radius: None,
        });
    }

    // Rectilinear staircases.
    for i in 0..count / 6 {
        let steps = 3 + (i % 6);
        let step = rng.range(0.5, 2.0);
        let mut pts = vec![(0.0, 0.0)];
        let (mut x, mut y) = (0.0, 0.0);
        for s in 0..steps {
            if s % 2 == 0 {
                y += step * rng.range(0.6, 2.0);
                pts.push((x, y));
                x += step * rng.range(0.6, 2.0);
                pts.push((x, y));
            } else {
                x += step * rng.range(0.6, 2.0);
                pts.push((x, y));
                y += step * rng.range(0.6, 2.0);
                pts.push((x, y));
            }
        }
        pts.push((x, 0.0));
        pts.push((0.0, 0.0));
        out.push(Case {
            name: format!("stair {i} steps={steps}"),
            kind: "synthetic",
            polygon: Polygon::new(ring(&pts), vec![]),
            exact_radius: None,
        });
    }

    // Square with 1..3 non-overlapping rectangular holes.
    for i in 0..count / 6 {
        let holes_n = 1 + (i % 3);
        let cell = 5.0;
        let mut holes = Vec::new();
        for h in 0..holes_n {
            let cx = 1.0 + cell * ((h % 3) as f64);
            let cy = 1.0 + cell * (((h / 3) % 3) as f64);
            let w = rng.range(1.0, 3.0);
            let hh = rng.range(1.0, 3.0);
            holes.push(ring(&[
                (cx, cy),
                (cx + w, cy),
                (cx + w, cy + hh),
                (cx, cy + hh),
                (cx, cy),
            ]));
        }
        out.push(Case {
            name: format!("holed {i} holes={holes_n}"),
            kind: "synthetic",
            polygon: Polygon::new(
                ring(&[
                    (0.0, 0.0),
                    (16.0, 0.0),
                    (16.0, 16.0),
                    (0.0, 16.0),
                    (0.0, 0.0),
                ]),
                holes,
            ),
            exact_radius: None,
        });
    }

    // Thin slivers, where the radius to diagonal ratio is brutal and pruning
    // has the least to work with.
    for i in 0..count / 6 {
        let w = rng.range(0.01, 0.4);
        let len = rng.range(5.0, 50.0);
        let skew = rng.range(0.0, 0.5) * len;
        out.push(Case {
            name: format!("sliver {i} w={w:.3}"),
            kind: "synthetic",
            polygon: Polygon::new(
                ring(&[
                    (0.0, 0.0),
                    (len, 0.0),
                    (len - skew, w),
                    (skew * 0.5, w * rng.range(0.2, 1.0)),
                    (0.0, 0.0),
                ]),
                vec![],
            ),
            exact_radius: None,
        });
    }

    // Random convex hulls.
    for i in 0..count / 6 {
        let n = 12 + (i % 40);
        let mut pts: Vec<(f64, f64)> = (0..n)
            .map(|_| (rng.range(-5.0, 5.0), rng.range(-5.0, 5.0)))
            .collect();
        pts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pts.dedup();
        let hull = convex_hull(&pts);
        if hull.len() >= 4 {
            out.push(Case {
                name: format!("hull {i} n={n}"),
                kind: "synthetic",
                polygon: Polygon::new(ring(&hull), vec![]),
                exact_radius: None,
            });
        }
    }

    out
}

fn real_cases(dirs: &[PathBuf]) -> (Vec<Case>, BTreeMap<String, usize>) {
    let mut cases = Vec::new();
    let mut skipped: BTreeMap<String, usize> = BTreeMap::new();
    let mut files: Vec<PathBuf> = Vec::new();
    for dir in dirs {
        let Ok(entries) = fs::read_dir(dir) else {
            *skipped
                .entry(format!("missing dir {}", dir.display()))
                .or_default() += 1;
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) == Some("wkt") {
                files.push(path);
            }
        }
    }
    files.sort();

    for path in files {
        let stem = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("?")
            .to_string();
        let Ok(text) = fs::read_to_string(&path) else {
            *skipped.entry(format!("{stem}: unreadable")).or_default() += 1;
            continue;
        };
        let parsed = match Wkt::<f64>::from_str(&text) {
            Ok(v) => v,
            Err(e) => {
                *skipped
                    .entry(format!("{stem}: parse error {e}"))
                    .or_default() += 1;
                continue;
            }
        };
        let geom = match geo_types::Geometry::<f64>::try_from(parsed) {
            Ok(g) => g,
            Err(e) => {
                *skipped
                    .entry(format!("{stem}: conversion error {e}"))
                    .or_default() += 1;
                continue;
            }
        };
        let components: Vec<Polygon<f64>> = match geom {
            geo_types::Geometry::Polygon(p) => vec![p],
            geo_types::Geometry::MultiPolygon(mp) => mp.0,
            other => {
                *skipped
                    .entry(format!("{stem}: non-polygonal ({other:?})"))
                    .or_default() += 1;
                continue;
            }
        };
        let multi = components.len() > 1;
        for (idx, poly) in components.into_iter().enumerate() {
            if cases.len() >= REAL_COMPONENT_LIMIT {
                *skipped
                    .entry("real component cap reached".to_string())
                    .or_default() += 1;
                continue;
            }
            let verts =
                poly.exterior().0.len() + poly.interiors().iter().map(|r| r.0.len()).sum::<usize>();
            if verts > COMPONENT_VERTEX_LIMIT {
                *skipped
                    .entry(format!("{stem}#{idx}: {verts} vertices over cap"))
                    .or_default() += 1;
                continue;
            }
            let name = if multi {
                format!("{stem}#{idx}")
            } else {
                stem.clone()
            };
            cases.push(Case {
                name,
                kind: "real",
                polygon: poly,
                exact_radius: None,
            });
        }
    }

    (cases, skipped)
}

struct Measurements {
    native: Option<MicResult>,
    geos: Option<MicResult>,
    native_err: Option<String>,
    geos_err: Option<String>,
    clearance_native: f64,
    clearance_geos: f64,
    us_native: f64,
    us_geos: f64,
    inside_native: bool,
    inside_geos: bool,
}

fn measure(poly: &Polygon<f64>, host: &HostPolygon) -> Measurements {
    let options = MicOptions::default();
    let tolerance = grid_tolerance(host);

    let native_start = Instant::now();
    let native = maximum_inscribed_circle(poly, &options);
    let us_native = native_start.elapsed().as_secs_f64() * 1e6;

    let geos_start = Instant::now();
    let geos = solve_mic_geos_with_tolerance(host, tolerance);
    let us_geos = geos_start.elapsed().as_secs_f64() * 1e6;

    let clearance = |r: Option<&MicResult>| match r {
        Some(res) => clearance_brute(poly, res.center.x(), res.center.y()),
        None => 0.0,
    };
    let inside = |r: Option<&MicResult>| match r {
        Some(res) => poly.contains(&Point::new(res.center.x(), res.center.y())),
        None => false,
    };

    let native_ok = native.as_ref().ok();
    let geos_ok = geos.as_ref().ok();
    let clearance_native = clearance(native_ok);
    let clearance_geos = clearance(geos_ok);
    let inside_native = inside(native_ok);
    let inside_geos = inside(geos_ok);
    let native_err = native.as_ref().err().map(|e| e.to_string());
    let geos_err = geos.as_ref().err().map(|e| e.to_string());

    Measurements {
        native: native.ok(),
        geos: geos.ok(),
        native_err,
        geos_err,
        clearance_native,
        clearance_geos,
        us_native,
        us_geos,
        inside_native,
        inside_geos,
    }
}

fn clearance_brute(poly: &Polygon<f64>, x: f64, y: f64) -> f64 {
    let mut best = f64::INFINITY;
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors().iter()) {
        for w in ring.0.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (dx, dy) = (b.x - a.x, b.y - a.y);
            let len2 = dx * dx + dy * dy;
            let t = if len2 <= 0.0 {
                0.0
            } else {
                (((x - a.x) * dx + (y - a.y) * dy) / len2).clamp(0.0, 1.0)
            };
            let px = a.x + t * dx;
            let py = a.y + t * dy;
            let d = ((x - px) * (x - px) + (y - py) * (y - py)).sqrt();
            if d < best {
                best = d;
            }
        }
    }
    if best.is_finite() {
        best
    } else {
        0.0
    }
}

fn sampling_floor(poly: &Polygon<f64>, host: &HostPolygon) -> Option<(f64, f64, f64)> {
    let verts = host.ring_coords(0).len()
        + (1..host.ring_count())
            .map(|i| host.ring_coords(i).len())
            .sum::<usize>();
    if verts > FLOOR_VERTEX_LIMIT {
        return None;
    }
    let (min_x, min_y, max_x, max_y) = host.bounds()?;
    let nb = NearestBoundaryIndex::new(SegmentIndex::from_host(host));
    let grid = if verts > 2_000 { 32 } else { FLOOR_GRID };
    let mut best = (0.0f64, min_x, min_y);
    for i in 1..grid {
        for j in 1..grid {
            let x = min_x + (max_x - min_x) * (i as f64) / (grid as f64);
            let y = min_y + (max_y - min_y) * (j as f64) / (grid as f64);
            if !poly.contains(&Point::new(x, y)) {
                continue;
            }
            if let Some((d, _)) = nb.nearest_distance_sq(x, y) {
                let d = d.max(0.0).sqrt();
                if d > best.0 {
                    best = (d, x, y);
                }
            }
        }
    }
    Some(best)
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dirs: Vec<PathBuf> = if args.is_empty() {
        default_wkt_dirs()
    } else {
        args.iter().map(PathBuf::from).collect()
    };

    let (real, skipped) = real_cases(&dirs);
    let mut cases = analytic_cases();
    cases.extend(real);
    cases.extend(synthetic_cases(150));

    println!(
        "corpus: {} cases (analytic {}, real {}, synthetic {})",
        cases.len(),
        cases.iter().filter(|c| c.kind == "analytic").count(),
        cases.iter().filter(|c| c.kind == "real").count(),
        cases.iter().filter(|c| c.kind == "synthetic").count()
    );
    if !skipped.is_empty() {
        println!("skipped inputs:");
        for (what, n) in &skipped {
            println!("  {n:>4}  {what}");
        }
    }

    let mut outcomes: Vec<Outcome> = Vec::new();
    for case in &cases {
        let host = match HostPolygon::from_polygon(&case.polygon) {
            Ok(h) => h,
            Err(e) => {
                outcomes.push(Outcome {
                    name: case.name.clone(),
                    kind: case.kind,
                    verts: case.polygon.exterior().0.len(),
                    native: None,
                    geos: None,
                    native_err: Some(format!("input rejected: {e}")),
                    geos_err: Some("input rejected".to_string()),
                    clearance_native: 0.0,
                    clearance_geos: 0.0,
                    inside_native: false,
                    inside_geos: false,
                    floor: None,
                    floor_proves_native: None,
                    floor_proves_geos: None,
                    engine_native: String::new(),
                    approach_native: String::new(),
                    engine_geos: String::new(),
                    candidates_native: 0,
                    candidates_geos: 0,
                    unclosed_rings: 0,
                    us_native: 0.0,
                    us_geos: 0.0,
                    exact_radius: case.exact_radius,
                });
                continue;
            }
        };

        let m = measure(&case.polygon, &host);

        let engine_native = m
            .native
            .as_ref()
            .map(|r| format!("{:?}", r.used_engine))
            .unwrap_or_default();
        let approach_native = m
            .native
            .as_ref()
            .map(|r| format!("{:?}", r.approach))
            .unwrap_or_default();
        let engine_geos = m
            .geos
            .as_ref()
            .map(|r| format!("{:?}", r.used_engine))
            .unwrap_or_default();
        let candidates_native = m.native.as_ref().map(|r| r.candidate_count).unwrap_or(0);
        let candidates_geos = m.geos.as_ref().map(|r| r.candidate_count).unwrap_or(0);

        let unclosed_rings = std::iter::once(case.polygon.exterior())
            .chain(case.polygon.interiors().iter())
            .filter(|ring| match (ring.0.first(), ring.0.last()) {
                (Some(a), Some(b)) => a.x != b.x || a.y != b.y,
                _ => true,
            })
            .count();

        let floor = sampling_floor(&case.polygon, &host);
        let (mut floor_proves_native, mut floor_proves_geos) = (None, None);
        if let Some((fd, fx, fy)) = floor {
            // Confirm the sampled point with plain arithmetic before it is
            // allowed to prove that a solver came back short.
            let confirmed = clearance_brute(&case.polygon, fx, fy);
            let fd = fd.min(confirmed).max(0.0);
            if let Some(r) = &m.native {
                if fd > r.radius * (1.0 + 1e-9) {
                    floor_proves_native = Some((fd - r.radius) / r.radius.max(1e-300));
                }
            }
            if let Some(r) = &m.geos {
                if fd > r.radius * (1.0 + 1e-9) {
                    floor_proves_geos = Some((fd - r.radius) / r.radius.max(1e-300));
                }
            }
        }

        outcomes.push(Outcome {
            name: case.name.clone(),
            kind: case.kind,
            verts: case.polygon.exterior().0.len()
                + case
                    .polygon
                    .interiors()
                    .iter()
                    .map(|r| r.0.len())
                    .sum::<usize>(),
            native: m.native,
            geos: m.geos,
            native_err: m.native_err,
            geos_err: m.geos_err,
            clearance_native: m.clearance_native,
            clearance_geos: m.clearance_geos,
            inside_native: m.inside_native,
            inside_geos: m.inside_geos,
            floor: floor.map(|(d, _, _)| d),
            floor_proves_native,
            floor_proves_geos,
            engine_native,
            approach_native,
            engine_geos,
            candidates_native,
            candidates_geos,
            unclosed_rings,
            us_native: m.us_native,
            us_geos: m.us_geos,
            exact_radius: case.exact_radius,
        });
    }

    // ---- aggregates ----
    let both: Vec<&Outcome> = outcomes
        .iter()
        .filter(|o| o.native.is_some() && o.geos.is_some())
        .collect();
    let native_errors = outcomes.iter().filter(|o| o.native_err.is_some()).count();
    let geos_errors = outcomes.iter().filter(|o| o.geos_err.is_some()).count();

    let mut gaps: Vec<(f64, String)> = Vec::new();
    for o in &both {
        let rn = o.native.as_ref().unwrap().radius;
        let rg = o.geos.as_ref().unwrap().radius;
        let denom = rn.abs().max(rg.abs()).max(1e-300);
        gaps.push(((rn - rg) / denom, o.name.clone()));
    }
    let mut abs_sorted: Vec<f64> = gaps.iter().map(|(g, _)| g.abs()).collect();
    abs_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let us_total_native: f64 = both.iter().map(|o| o.us_native).sum();
    let us_total_geos: f64 = both.iter().map(|o| o.us_geos).sum();

    println!(
        "\n=== radius parity, native vs GEOS ({} cases solved by both) ===",
        both.len()
    );
    println!("cases where the input was rejected: native {native_errors}  geos {geos_errors}");
    if !abs_sorted.is_empty() {
        println!(
            "relative gap   max {:.3e}   p90 {:.3e}   median {:.3e}   mean {:.3e}",
            abs_sorted[abs_sorted.len() - 1],
            percentile(&abs_sorted, 0.9),
            percentile(&abs_sorted, 0.5),
            abs_sorted.iter().sum::<f64>() / abs_sorted.len() as f64,
        );
        for t in [1e-15, 1e-12, 1e-9, 1e-7, 1e-6, 1e-4] {
            let n = abs_sorted.iter().filter(|g| **g > t).count();
            println!("  gaps above {t:.0e}: {n}");
        }
    }

    let mut native_smaller: Vec<(f64, String)> =
        gaps.iter().filter(|(g, _)| *g < -1e-9).cloned().collect();
    native_smaller.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut geos_smaller: Vec<(f64, String)> =
        gaps.iter().filter(|(g, _)| *g > 1e-9).cloned().collect();
    geos_smaller.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    println!(
        "\n  native returned the smaller circle on {}   GEOS on {}",
        native_smaller.len(),
        geos_smaller.len()
    );
    println!("  worst cases where native came back smaller:");
    for (g, name) in native_smaller.iter().take(8) {
        println!("    {g:+.3e}  {name}");
    }
    println!("  worst cases where GEOS came back smaller:");
    for (g, name) in geos_smaller.iter().take(8) {
        println!("    {g:+.3e}  {name}");
    }

    // ---- independent checks ----
    println!("\n=== independent checks (no solver involved) ===");
    let mut over_native = Vec::new();
    let mut over_geos = Vec::new();
    let mut outside_native = Vec::new();
    let mut outside_geos = Vec::new();
    let mut floor_native = Vec::new();
    let mut floor_geos = Vec::new();
    for o in &outcomes {
        if let Some(r) = &o.native {
            let over = r.radius - o.clearance_native;
            if over > 1e-9 * r.radius.max(1.0) {
                over_native.push((
                    over / r.radius.max(1e-300),
                    format!("{} [{}]", o.name, o.engine_native),
                ));
            }
            if !o.inside_native {
                outside_native.push(o.name.clone());
            }
        }
        if let Some(r) = &o.geos {
            let over = r.radius - o.clearance_geos;
            if over > 1e-9 * r.radius.max(1.0) {
                over_geos.push((
                    over / r.radius.max(1e-300),
                    format!("{} [{}]", o.name, o.engine_geos),
                ));
            }
            if !o.inside_geos {
                outside_geos.push(o.name.clone());
            }
        }
        if let Some(def) = o.floor_proves_native {
            floor_native.push((def, format!("{} [{}]", o.name, o.engine_native)));
        }
        if let Some(def) = o.floor_proves_geos {
            floor_geos.push((def, format!("{} [{}]", o.name, o.engine_geos)));
        }
    }
    over_native.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    over_geos.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    floor_native.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    floor_geos.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    println!(
        "radius above the true clearance (circle does not fit): native {}   GEOS {}",
        over_native.len(),
        over_geos.len()
    );
    for (over, name) in over_native.iter().take(5) {
        println!("    native over by {over:+.3e}  {name}");
    }
    for (over, name) in over_geos.iter().take(5) {
        println!("    GEOS   over by {over:+.3e}  {name}");
    }
    println!(
        "returned center outside the polygon: native {}   GEOS {}",
        outside_native.len(),
        outside_geos.len()
    );
    for name in outside_native.iter().take(5) {
        println!("    native {name}");
    }
    for name in outside_geos.iter().take(5) {
        println!("    GEOS   {name}");
    }
    let unclosed_cases = outcomes.iter().filter(|o| o.unclosed_rings > 0).count();
    println!(
        "cases with an unclosed ring (the check would then miss one segment): {unclosed_cases}"
    );
    for o in outcomes.iter().filter(|o| o.unclosed_rings > 0).take(5) {
        println!("    {} ({} rings)", o.name, o.unclosed_rings);
    }
    let with_floor = outcomes.iter().filter(|o| o.floor.is_some()).count();
    println!(
        "sampled floor above a returned radius ({with_floor} cases sampled): native {}   GEOS {}",
        floor_native.len(),
        floor_geos.len()
    );
    for (def, name) in floor_native.iter().take(8) {
        let cands = outcomes
            .iter()
            .find(|o| name.starts_with(&o.name))
            .map(|o| o.candidates_native)
            .unwrap_or(0);
        println!("    native short by {def:.3e}  {name}  candidates={cands}");
    }
    for (def, name) in floor_geos.iter().take(8) {
        println!("    GEOS   short by {def:.3e}  {name}");
    }

    // ---- closed-form check ----
    println!("\n=== closed-form shapes ===");
    println!(
        "{:<24} {:>14} {:>13} {:>13}",
        "shape", "exact radius", "native rel", "GEOS rel"
    );
    for o in outcomes.iter().filter(|o| o.exact_radius.is_some()) {
        let exact = o.exact_radius.unwrap();
        let rel = |r: Option<f64>| match r {
            Some(v) => format!("{:+.3e}", (v - exact) / exact),
            None => "error".to_string(),
        };
        println!(
            "{:<24} {:>14.9} {:>13} {:>13}",
            o.name,
            exact,
            rel(o.native.as_ref().map(|r| r.radius)),
            rel(o.geos.as_ref().map(|r| r.radius))
        );
    }

    // ---- timing ----
    println!("\n=== wall time over the shared corpus ===");
    println!(
        "native total {:.3}s   GEOS total {:.3}s   ratio {:.2}x",
        us_total_native / 1e6,
        us_total_geos / 1e6,
        us_total_geos / us_total_native.max(1e-9)
    );
    for kind in ["analytic", "real", "synthetic"] {
        let sel: Vec<&&Outcome> = both.iter().filter(|o| o.kind == kind).collect();
        if sel.is_empty() {
            continue;
        }
        let n: f64 = sel.iter().map(|o| o.us_native).sum();
        let g: f64 = sel.iter().map(|o| o.us_geos).sum();
        println!(
            "  {:<10} {:>4} cases   native {:>9.3}s   GEOS {:>9.3}s   {:.2}x",
            kind,
            sel.len(),
            n / 1e6,
            g / 1e6,
            g / n.max(1e-9)
        );
    }

    // ---- csv ----
    let out_path = PathBuf::from("target/parity_geos_report.csv");
    let mut csv = String::from(
        "name,kind,verts,radius_native,radius_geos,rel_gap,clearance_native,clearance_geos,floor,inside_native,inside_geos,us_native,us_geos,native_err,geos_err,engine_native,approach_native,engine_geos,candidates_native,candidates_geos,radius_exact\n",
    );
    for o in &outcomes {
        let rn = o.native.as_ref().map(|r| r.radius);
        let rg = o.geos.as_ref().map(|r| r.radius);
        let gap = match (rn, rg) {
            (Some(a), Some(b)) => (a - b) / a.abs().max(b.abs()).max(1e-300),
            _ => f64::NAN,
        };
        csv.push_str(&format!(
            "{},{},{},{},{},{:.6e},{:.9e},{:.9e},{},{},{},{:.1},{:.1},{},{},{},{},{},{},{},{}\n",
            o.name.replace(',', ";"),
            o.kind,
            o.verts,
            rn.map(|v| format!("{v:.9}")).unwrap_or_default(),
            rg.map(|v| format!("{v:.9}")).unwrap_or_default(),
            gap,
            o.clearance_native,
            o.clearance_geos,
            o.floor.map(|f| format!("{f:.9}")).unwrap_or_default(),
            o.inside_native,
            o.inside_geos,
            o.us_native,
            o.us_geos,
            o.native_err.as_deref().unwrap_or(""),
            o.geos_err.as_deref().unwrap_or(""),
            o.engine_native,
            o.approach_native,
            o.engine_geos,
            o.candidates_native,
            o.candidates_geos,
            o.exact_radius
                .map(|r| format!("{r:.17e}"))
                .unwrap_or_default(),
        ));
    }
    if let Some(parent) = out_path.parent() {
        let _ = fs::create_dir_all(parent);
    }
    match fs::write(&out_path, csv) {
        Ok(()) => println!("\nper-case csv: {}", out_path.display()),
        Err(e) => println!("\ncsv write failed: {e}"),
    }
}
