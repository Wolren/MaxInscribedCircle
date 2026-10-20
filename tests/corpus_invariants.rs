//! Corpus-level invariants that run without GEOS.
//!
//! The GEOS parity harness lives in `examples/parity_geos.rs` because it links
//! the C library. The checks here need no external solver: every returned
//! circle is re-measured against the segment table, and a dense sample of
//! interior points proves the answer is not below the floor. Both properties
//! hold for any correct solver, so they gate every push.

use geo_types::{Coord, LineString, Polygon};
use max_inscribed_circle::mic::input::{HostPolygon, SegmentIndex};
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

fn poly(pts: &[(f64, f64)]) -> Polygon<f64> {
    let mut v: Vec<Coord<f64>> = pts.iter().map(|p| Coord { x: p.0, y: p.1 }).collect();
    v.push(v[0]);
    Polygon::new(LineString::from(v), vec![])
}

fn clearance(p: &Polygon<f64>, x: f64, y: f64) -> f64 {
    let host = HostPolygon::from_polygon(p).expect("host");
    let seg = SegmentIndex::from_host(&host);
    seg.batch_point_segment_distance_sq(x, y, 0, seg.len())
        .sqrt()
}

/// Lower bound on the MIC: the best clearance over a grid of interior points.
fn floor(p: &Polygon<f64>, steps: usize) -> f64 {
    let host = HostPolygon::from_polygon(p).expect("host");
    let pip = max_inscribed_circle::mic::index::PipIndex::new(&host);
    let seg = SegmentIndex::from_host(&host);
    let (min_x, min_y, max_x, max_y) = host.bounds().expect("bounds");
    let mut best = 0.0f64;
    for i in 0..=steps {
        let x = min_x + (max_x - min_x) * (i as f64) / (steps as f64);
        for j in 0..=steps {
            let y = min_y + (max_y - min_y) * (j as f64) / (steps as f64);
            if !pip.contains_strict_xy(x, y) {
                continue;
            }
            let d = seg
                .batch_point_segment_distance_sq(x, y, 0, seg.len())
                .sqrt();
            if d > best {
                best = d;
            }
        }
    }
    best
}

struct Case {
    name: &'static str,
    polygon: Polygon<f64>,
}

fn corpus() -> Vec<Case> {
    let mut cases = Vec::new();

    let closed_form: Vec<(&'static str, Polygon<f64>)> = vec![
        (
            "square 10",
            poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        ),
        (
            "rect 8x2",
            poly(&[(0.0, 0.0), (8.0, 0.0), (8.0, 2.0), (0.0, 2.0)]),
        ),
        (
            "tri equilateral",
            poly(&[(0.0, 0.0), (10.0, 0.0), (5.0, 8.660_254_038)]),
        ),
        (
            "L shape",
            poly(&[
                (0.0, 0.0),
                (4.0, 0.0),
                (4.0, 1.0),
                (1.0, 1.0),
                (1.0, 3.0),
                (0.0, 3.0),
            ]),
        ),
        (
            "U shape",
            poly(&[
                (0.0, 0.0),
                (10.0, 0.0),
                (10.0, 1.0),
                (6.0, 1.0),
                (6.0, 9.0),
                (4.0, 9.0),
                (4.0, 1.0),
                (0.0, 1.0),
            ]),
        ),
        (
            "thin corridor",
            poly(&[(0.0, 0.0), (100.0, 0.0), (100.0, 2.0), (0.0, 2.0)]),
        ),
        (
            "notched corridor",
            poly(&[
                (0.0, 0.0),
                (100.0, 0.0),
                (100.0, 2.0),
                (51.0, 2.0),
                (51.0, 3.0),
                (50.0, 3.0),
                (50.0, 2.0),
                (0.0, 2.0),
            ]),
        ),
        (
            "staircase",
            poly(&[
                (0.0, 0.0),
                (6.0, 0.0),
                (6.0, 1.0),
                (4.0, 1.0),
                (4.0, 2.0),
                (2.0, 2.0),
                (2.0, 3.0),
                (0.0, 3.0),
            ]),
        ),
        (
            "concave star-ish",
            poly(&[
                (0.0, 0.0),
                (4.0, -1.0),
                (8.0, 0.0),
                (7.0, 3.0),
                (4.0, 2.0),
                (1.0, 3.0),
            ]),
        ),
    ];
    for (name, polygon) in closed_form {
        cases.push(Case { name, polygon });
    }

    // Regular n-gons: MIC is the inradius, so the answer is known exactly.
    for n in [5usize, 12, 64, 256] {
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
                (10.0 * a.cos(), 10.0 * a.sin())
            })
            .collect();
        let name: &'static str = match n {
            5 => "regular 5-gon",
            12 => "regular 12-gon",
            64 => "regular 64-gon",
            _ => "regular 256-gon",
        };
        cases.push(Case {
            name,
            polygon: poly(&pts),
        });
    }

    cases
}

#[test]
fn every_returned_circle_fits_its_polygon() {
    let opts = MicOptions::default();
    for case in corpus() {
        let r = maximum_inscribed_circle(&case.polygon, &opts)
            .unwrap_or_else(|e| panic!("{}: solve failed: {e}", case.name));
        let measured = clearance(&case.polygon, r.center.x(), r.center.y());
        let rel = (r.radius - measured).abs() / measured.max(f64::MIN_POSITIVE);
        assert!(
            rel <= 1e-9,
            "{}: reported {} but the segment table measures {}",
            case.name,
            r.radius,
            measured
        );
    }
}

#[test]
fn no_answer_falls_below_a_sampled_floor() {
    let opts = MicOptions::default();
    for case in corpus() {
        let r = maximum_inscribed_circle(&case.polygon, &opts).expect("solve");
        let sampled = floor(&case.polygon, 96);
        let slack = 1e-6 * sampled.max(1.0);
        assert!(
            r.radius + slack >= sampled,
            "{}: returned {} but an interior point already clears {}",
            case.name,
            r.radius,
            sampled
        );
    }
}

#[test]
fn regular_polygons_match_the_inradius() {
    let opts = MicOptions::default();
    for n in [5usize, 12, 64, 256] {
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
                (10.0 * a.cos(), 10.0 * a.sin())
            })
            .collect();
        let p = poly(&pts);
        let r = maximum_inscribed_circle(&p, &opts).expect("solve");
        let truth = 10.0 * (std::f64::consts::PI / n as f64).cos();
        assert!(
            (r.radius - truth).abs() / truth < 1e-6,
            "{n}-gon: got {} want {truth}",
            r.radius
        );
    }
}
