//! Maximum Inscribed Circle (MIC) solvers for polygonal inputs.
//!
//! One pipeline: closed-form fast paths for the shape classes where the
//! optimum is known, then the grid solver for everything else. There is no
//! external solver in the pipeline; the GEOS comparison lives behind the
//! `bench-geos` feature in `solver::bench_geos` and is never called from here.

pub mod index;
pub mod input;
pub mod solver;
pub mod workspace;

use crate::tuning;
use geo_types::{MultiPolygon, Point, Polygon};
use rayon::prelude::*;
use thiserror::Error;

use self::input::HostPolygon;
use self::solver::grid::solve_grid;
use self::workspace::MicWorkspace;

/// Solver configuration.
///
/// The solve has no tunable knobs: the grid tolerance scales with the input
/// diagonal and the fast paths are closed form. The type is kept so callers
/// pass options by reference and future knobs stay source compatible.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MicOptions {}

/// Solver that produced the final result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MicUsedEngine {
    /// The grid solver produced the circle.
    Grid,
    /// The GEOS comparison solver produced the circle, behind `bench-geos`.
    Geos,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Solver path that produced the result.
pub enum MicApproach {
    /// Analytical incenter fast path for a triangular ring.
    FastTriangle,
    /// Fast path for an axis-aligned rectangular ring.
    FastRectangle,
    /// Chebyshev-center fast path for a convex quadrilateral ring.
    FastConvexQuad,
    /// Fast path for a convex n-gon, solved by half-plane feasibility search.
    FastConvexN,
    /// Uniform grid search over interior cell centers.
    GridSolver,
    /// GEOS maximum inscribed circle, behind the `bench-geos` feature.
    Geos,
}

/// MIC solve result.
#[derive(Debug, Clone)]
pub struct MicResult {
    /// Circle center point.
    pub center: Point<f64>,
    /// Circle radius in input coordinate units.
    pub radius: f64,
    /// Squared radius, radius * radius.
    pub radius_sq: f64,
    /// Indices of the boundary segments the circle touches.
    pub support_segments: Vec<usize>,
    /// Number of candidate centers evaluated during the solve.
    pub candidate_count: usize,
    /// Solver that produced this result.
    pub used_engine: MicUsedEngine,
    /// Index of the multipolygon component that holds the circle, if any.
    pub component_index: Option<usize>,
    /// Solver path that produced this result.
    pub approach: MicApproach,
}

/// MIC solver error.
#[derive(Debug, Error)]
pub enum MicError {
    /// Input polygon is malformed, empty, or degenerate.
    #[error("invalid polygon input: {0}")]
    InvalidInput(String),
    /// No interior point was found, so no circle exists.
    #[error("no valid MIC candidate found")]
    NoCircleFound,
    /// The GEOS comparison solve failed, behind the `bench-geos` feature.
    #[error("GEOS comparison solve failed: {0}")]
    GeosFailed(String),
    /// The GEOS result could not be parsed into a circle.
    #[error("unsupported GEOS MIC output: {0}")]
    UnsupportedGeosOutput(String),
}

/// Solve MIC on a single polygon.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn maximum_inscribed_circle(
    poly: &Polygon<f64>,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    let host = HostPolygon::from_polygon(poly)?;
    solve_on_host_polygon(&host, opts)
}

/// Solve MIC reusing a workspace, which holds the indexes built for its polygon.
///
/// The workspace owns the polygon it was built for, so this solves that
/// polygon. Build a fresh [`MicWorkspace`] to solve a different one.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn maximum_inscribed_circle_with_workspace(
    workspace: &mut MicWorkspace,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    solve_with_indexes(
        &workspace.host,
        &workspace.pip_index,
        &workspace.nb_index,
        opts,
    )
}

/// Solve MIC on a multipolygon by solving each component and keeping the best.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn maximum_inscribed_circle_multipolygon(
    multi: &MultiPolygon<f64>,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    if multi.0.is_empty() {
        return Err(MicError::InvalidInput(
            "multipolygon has no components".to_string(),
        ));
    }

    let mut best: Option<MicResult> = None;
    let mut last_error: Option<MicError> = None;

    for (idx, poly) in multi.0.iter().enumerate() {
        match maximum_inscribed_circle(poly, opts) {
            Ok(mut result) => {
                result.component_index = Some(idx);
                let replace = best
                    .as_ref()
                    .map(|current| result.radius_sq > current.radius_sq)
                    .unwrap_or(true);
                if replace {
                    best = Some(result);
                }
            }
            Err(err) => {
                last_error = Some(err);
            }
        }
    }

    best.ok_or_else(|| last_error.unwrap_or(MicError::NoCircleFound))
}

/// Solve MIC on many polygons in parallel using rayon.
///
/// Each polygon is solved independently on a rayon thread pool.
/// Returns results in the same order as input.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn maximum_inscribed_circle_slice(
    polys: &[Polygon<f64>],
    opts: &MicOptions,
) -> Vec<Result<MicResult, MicError>> {
    polys
        .par_iter()
        .map(|poly| maximum_inscribed_circle(poly, opts))
        .collect()
}

fn solve_on_host_polygon(host: &HostPolygon, _opts: &MicOptions) -> Result<MicResult, MicError> {
    // Phase 0: Try analytical fast path for simple shapes (triangle, rectangle, convex quad).
    // We intentionally do NOT dispatch arbitrary N-gons to fast_convex_n here:
    // real-world data has near-collinear reflex vertices and numerical noise that
    // can fool convexity checks, producing kernel circles smaller than the true MIC.
    if host.ring_count() == 1 {
        let outer_len = host.outer_ring().len();
        let result = if outer_len == 4 {
            // Triangle - exact incenter formula
            solver::fast::fast_triangle(host)
        } else if outer_len == 5 {
            // Rectangle (cheapest O(1)), then general convex quad (LP solver)
            solver::fast::fast_rectangle(host).or_else(|| solver::fast::fast_convex_quad(host))
        } else {
            None
        };
        if let Some(mut result) = result {
            // Every radius this function returns is the distance the segment
            // table reports at the returned centre, so the reported circle
            // always fits. The fast paths solve the shape class they recognise;
            // this measurement confirms the polygon really is that class,
            // because near-collinear reflex vertices in real data can fool the
            // convexity tests. A shape that only looks like the class falls
            // through to the grid instead of returning a radius it cannot
            // support.
            let seg_idx = input::SegmentIndex::from_host(host);
            let actual_sq = seg_idx.batch_point_segment_distance_sq(
                result.center.x(),
                result.center.y(),
                0,
                seg_idx.len(),
            );
            let actual = actual_sq.max(0.0).sqrt();
            let relative = if actual > 0.0 {
                (result.radius - actual).abs() / actual
            } else {
                f64::INFINITY
            };
            if relative <= tuning::MIC_FAST_PATH_REL_TOL {
                result.radius = actual;
                result.radius_sq = actual_sq;
                return Ok(result);
            }
        }
    }

    // Phase 1: grid search over the bounding box, which answers every shape
    // the fast paths did not claim.
    let workspace = MicWorkspace::new(host.clone())?;
    let tolerance = grid_tolerance(host);

    if let Some(grid_result) =
        solve_grid(host, &workspace.pip_index, &workspace.nb_index, tolerance)
    {
        return Ok(grid_result);
    }

    // The seeds cover the whole box, so reaching here means the polygon has no
    // interior point at all, which a zero-area ring catches earlier. Report it
    // rather than return a circle the input does not contain.
    Err(MicError::NoCircleFound)
}

/// Solve with indexes the caller already built, skipping index construction.
fn solve_with_indexes(
    host: &HostPolygon,
    pip: &index::PipIndex,
    nb: &index::NearestBoundaryIndex,
    _opts: &MicOptions,
) -> Result<MicResult, MicError> {
    let tolerance = grid_tolerance(host);
    solve_grid(host, pip, nb, tolerance).ok_or(MicError::NoCircleFound)
}

/// Absolute refinement tolerance for the grid stage on a host polygon:
/// `diag * 1e-7`, where `diag` is the bounding box diagonal as it stands.
///
/// Exposed so a comparison harness can hand another solver the same tolerance
/// the native pipeline uses, instead of comparing solvers at mixed tolerances.
///
/// The diagonal carries no floor, because a floor would tie the accuracy of the
/// answer to the coordinate units of the input: a polygon a few hundred metres
/// wide stored in degrees has a diagonal near `1e-3`, so a floor of one unit
/// makes the tolerance a whole percent of its own inscribed radius and the
/// answer stops meaning anything. Scale invariance here means the same polygon
/// in degrees and in metres converges to the same circle.
pub fn grid_tolerance(host: &HostPolygon) -> f64 {
    let Some((min_x, min_y, max_x, max_y)) = host.bounds() else {
        return 0.0;
    };
    let diag = (max_x - min_x).hypot(max_y - min_y);
    diag * 1e-7
}

#[cfg(test)]
mod tests {
    use super::*;
    use geo_types::{Coord, LineString, Polygon};

    fn square(size: f64) -> Polygon<f64> {
        Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord { x: size, y: 0.0 },
                Coord { x: size, y: size },
                Coord { x: 0.0, y: size },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        )
    }

    #[test]
    fn mic_exact_square_radius() {
        let poly = square(10.0);
        let opts = MicOptions::default();
        let r = maximum_inscribed_circle(&poly, &opts).expect("solve failed");
        assert!((r.radius - 5.0).abs() < 1e-6, "radius={}", r.radius);
        assert!((r.center.x() - 5.0).abs() < 1e-6);
        assert!((r.center.y() - 5.0).abs() < 1e-6);
    }

    #[test]
    fn mic_exact_rectangle_radius_equal_half_width() {
        let poly = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord { x: 8.0, y: 0.0 },
                Coord { x: 8.0, y: 2.0 },
                Coord { x: 0.0, y: 2.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        );
        let opts = MicOptions::default();
        let r = maximum_inscribed_circle(&poly, &opts).expect("solve failed");
        assert!((r.radius - 1.0).abs() < 1e-6, "radius={}", r.radius);
    }

    #[test]
    fn mic_exact_triangle_center_inside() {
        let poly = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord { x: 10.0, y: 0.0 },
                Coord { x: 5.0, y: 10.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        );
        let opts = MicOptions::default();
        let r = maximum_inscribed_circle(&poly, &opts).expect("solve failed");
        assert!(r.radius > 1.0, "radius={}", r.radius);
        use geo::Contains;
        assert!(poly.contains(&r.center));
    }

    #[test]
    fn mic_degenerate_zero_area_returns_err() {
        let poly = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord { x: 1.0, y: 0.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        );
        let opts = MicOptions::default();
        assert!(maximum_inscribed_circle(&poly, &opts).is_err());
    }

    fn poly(pts: &[(f64, f64)]) -> Polygon<f64> {
        let mut v: Vec<Coord> = pts.iter().map(|p| Coord { x: p.0, y: p.1 }).collect();
        v.push(v[0]);
        Polygon::new(LineString::from(v), vec![])
    }

    fn assert_fits(poly: &Polygon<f64>, r: &MicResult) {
        use crate::mic::input::{HostPolygon, SegmentIndex};
        let host = HostPolygon::from_polygon(poly).expect("host");
        let seg = SegmentIndex::from_host(&host);
        let measured = seg
            .batch_point_segment_distance_sq(r.center.x(), r.center.y(), 0, seg.len())
            .sqrt();
        let rel = (r.radius - measured).abs() / measured.max(1e-300);
        assert!(
            rel <= 1e-9,
            "reported radius {} differs from measured clearance {}",
            r.radius,
            measured
        );
    }

    #[test]
    fn mic_rejects_nan_coordinate() {
        let bad = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord {
                    x: f64::NAN,
                    y: 0.0,
                },
                Coord { x: 10.0, y: 10.0 },
                Coord { x: 0.0, y: 10.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        );
        let err = maximum_inscribed_circle(&bad, &MicOptions::default()).unwrap_err();
        assert!(err.to_string().contains("NaN"), "got {err}");
    }

    #[test]
    fn mic_rejects_infinite_coordinate() {
        let bad = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord {
                    x: f64::INFINITY,
                    y: 0.0,
                },
                Coord { x: 10.0, y: 10.0 },
                Coord { x: 0.0, y: 10.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![],
        );
        assert!(maximum_inscribed_circle(&bad, &MicOptions::default()).is_err());
    }

    #[test]
    fn mic_respects_hole() {
        // 10x10 frame with a 4x4 hole: the circle must clear the hole, not
        // span the frame. Analytic value is sqrt(2) * (5 - 2) / 1 .. = 1.75736.
        let frame = Polygon::new(
            LineString::from(vec![
                Coord { x: 0.0, y: 0.0 },
                Coord { x: 10.0, y: 0.0 },
                Coord { x: 10.0, y: 10.0 },
                Coord { x: 0.0, y: 10.0 },
                Coord { x: 0.0, y: 0.0 },
            ]),
            vec![LineString::from(vec![
                Coord { x: 3.0, y: 3.0 },
                Coord { x: 7.0, y: 3.0 },
                Coord { x: 7.0, y: 7.0 },
                Coord { x: 3.0, y: 7.0 },
                Coord { x: 3.0, y: 3.0 },
            ])],
        );
        let r = maximum_inscribed_circle(&frame, &MicOptions::default()).expect("solve");
        assert!(
            (r.radius - 1.757_359_312_8).abs() < 1e-6,
            "radius={} expected the hole-limited value",
            r.radius
        );
        assert_fits(&frame, &r);
    }

    #[test]
    fn mic_concave_shape_fits() {
        let l = poly(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ]);
        let r = maximum_inscribed_circle(&l, &MicOptions::default()).expect("solve");
        assert!(
            (r.radius - 0.585_786_437_6).abs() < 1e-6,
            "radius={}",
            r.radius
        );
        assert_fits(&l, &r);
    }

    #[test]
    fn mic_thin_corridor_is_fast_and_fits() {
        // A ridge of equal clearance used to cost millions of cells here.
        let corridor = poly(&[(0.0, 0.0), (1000.0, 0.0), (1000.0, 2.0), (0.0, 2.0)]);
        let started = std::time::Instant::now();
        let r = maximum_inscribed_circle(&corridor, &MicOptions::default()).expect("solve");
        let elapsed = started.elapsed();
        assert!((r.radius - 1.0).abs() < 1e-9, "radius={}", r.radius);
        assert_fits(&corridor, &r);
        assert!(
            elapsed < std::time::Duration::from_millis(500),
            "thin corridor took {elapsed:?}"
        );
    }

    #[test]
    fn mic_notched_corridor_fits() {
        // The optimum sits in the notch basin, which the candidate stage alone
        // does not reach. Whatever stage answers, the circle must fit.
        let notched = poly(&[
            (0.0, 0.0),
            (100.0, 0.0),
            (100.0, 2.0),
            (51.0, 2.0),
            (51.0, 3.0),
            (50.0, 3.0),
            (50.0, 2.0),
            (0.0, 2.0),
        ]);
        let r = maximum_inscribed_circle(&notched, &MicOptions::default()).expect("solve");
        assert!(r.radius > 1.0, "radius={} lost the notch", r.radius);
        assert_fits(&notched, &r);
    }

    #[test]
    fn mic_multipolygon_reports_component() {
        let small = poly(&[(20.0, 20.0), (22.0, 20.0), (22.0, 22.0), (20.0, 22.0)]);
        let big = poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]);
        let multi = MultiPolygon(vec![small, big]);
        let r =
            maximum_inscribed_circle_multipolygon(&multi, &MicOptions::default()).expect("solve");
        assert!((r.radius - 5.0).abs() < 1e-6);
        assert_eq!(r.component_index, Some(1));
    }

    #[test]
    fn mic_batch_preserves_order() {
        let shapes: Vec<Polygon<f64>> = (1..=4)
            .map(|k| {
                poly(&[
                    (0.0, 0.0),
                    (k as f64, 0.0),
                    (k as f64, k as f64),
                    (0.0, k as f64),
                ])
            })
            .collect();
        let out = maximum_inscribed_circle_slice(&shapes, &MicOptions::default());
        assert_eq!(out.len(), 4);
        for (k, r) in out.iter().enumerate() {
            let side = (k + 1) as f64;
            assert!((r.as_ref().expect("solve").radius - side / 2.0).abs() < 1e-6);
        }
    }
}
