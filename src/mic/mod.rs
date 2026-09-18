//! Maximum Inscribed Circle (MIC) solvers for polygonal inputs.
//!
//! One pipeline: analytical fast paths, then the grid search, then the exact
//! candidate enumeration. There is no external solver in the pipeline; the
//! GEOS comparison lives behind the `bench-geos` feature in
//! `solver::bench_geos` and is never called from here.

pub mod index;
pub mod input;
pub mod solver;
pub mod workspace;

use geo_types::{MultiPolygon, Point, Polygon};
use rayon::prelude::*;
use thiserror::Error;

use self::input::HostPolygon;
use self::solver::exact::solve_exact;
use self::solver::grid::solve_grid;
use self::workspace::MicWorkspace;

/// Numeric robustness mode for the exact engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RobustMode {
    /// Fast finite-precision mode.
    FastF64,
    /// Extra candidate filtering and certification.
    Filtered,
}

/// Solver configuration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MicOptions {
    /// Robustness mode: Filtered adds candidate filtering to the exact solver.
    pub robust_mode: RobustMode,
    /// Use a BVH instead of a uniform grid for nearest-boundary queries.
    pub use_bvh: bool,
    /// Reserved and unused; batch parallel solves are provided by
    /// [`maximum_inscribed_circle_slice`].
    pub use_parallel: bool,
}

impl Default for MicOptions {
    fn default() -> Self {
        Self {
            robust_mode: RobustMode::Filtered,
            use_bvh: false,
            use_parallel: false,
        }
    }
}

/// Solver that produced the final result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MicUsedEngine {
    /// The native exact solver produced the circle.
    Exact,
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
    /// Exact candidate enumeration over segment and vertex Voronoi sites.
    VoronoiExact,
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
    /// No valid candidate center was found.
    #[error("no valid MIC candidate found")]
    NoCircleFound,
    /// The exact solver failed.
    #[error("exact MIC solver failed: {0}")]
    ExactFailed(String),
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

/// Solve MIC with a reusable workspace (avoids rebuilding indexes per call).
///
/// The workspace is rebuilt only if `host` changes; for repeated calls on
/// different polygons, create a fresh [`MicWorkspace`] each time.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn maximum_inscribed_circle_with_workspace(
    workspace: &mut MicWorkspace,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    solve_exact(workspace, opts).map_err(|err| MicError::ExactFailed(err.to_string()))
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

fn solve_on_host_polygon(host: &HostPolygon, opts: &MicOptions) -> Result<MicResult, MicError> {
    // Phase 0: Try analytical fast path for simple shapes (triangle, rectangle, convex quad).
    // We intentionally do NOT dispatch arbitrary N-gons to fast_convex_n here:
    // real-world data has near-collinear reflex vertices and numerical noise that
    // can fool convexity checks, producing kernel circles smaller than the true MIC.
    if host.ring_count() == 1 {
        let outer_len = host.outer_ring().len();
        let result = if outer_len == 4 {
            // Triangle - exact incenter formula
            solver::exact::fast_triangle(host)
        } else if outer_len == 5 {
            // Rectangle (cheapest O(1)), then general convex quad (LP solver)
            solver::exact::fast_rectangle(host).or_else(|| solver::exact::fast_convex_quad(host))
        } else {
            None
        };
        if let Some(result) = result {
            // Triangle and rectangle are exact - no verification needed.
            if outer_len <= 4 {
                return Ok(result);
            }
            // Convex quad: verify because the convexity check can be fooled by
            // nearly-collinear edges in real-world data.
            let seg_idx = input::SegmentIndex::from_host(host);
            let mut actual_sq = f64::INFINITY;
            for idx in 0..seg_idx.len() {
                let d =
                    seg_idx.point_segment_distance_sq(idx, result.center.x(), result.center.y());
                if d < actual_sq {
                    actual_sq = d;
                }
            }
            let actual = actual_sq.sqrt();
            if actual > 0.0 && (result.radius - actual).abs() / actual < 0.01 {
                return Ok(result);
            }
        }
    }

    // Phase 1: grid search over the bounding box. It answers most complex
    // shapes and returns before the exact stage runs.
    let workspace = MicWorkspace::new(host.clone(), opts.use_bvh)?;
    let tolerance = grid_tolerance(host);

    if let Some(grid_result) =
        solve_grid(host, &workspace.pip_index, &workspace.nb_index, tolerance)
    {
        return Ok(grid_result);
    }

    // Phase 2: exact candidate enumeration over Voronoi sites.
    run_exact(host, opts)
}

/// Tolerance the grid stage refines to on a host polygon: `max(diag * 1e-7, 1e-12)`.
///
/// Exposed so a comparison harness can hand another solver the same tolerance
/// the native pipeline uses, instead of comparing solvers at mixed tolerances.
pub fn grid_tolerance(host: &HostPolygon) -> f64 {
    let bounds = host.bounds().unwrap_or((0.0, 0.0, 1.0, 1.0));
    let diag = (bounds.2 - bounds.0).hypot(bounds.3 - bounds.1).max(1.0);
    (diag * 1e-7).max(1e-12)
}

fn run_exact(host: &HostPolygon, opts: &MicOptions) -> Result<MicResult, MicError> {
    let mut workspace = MicWorkspace::new(host.clone(), opts.use_bvh)?;
    solve_exact(&mut workspace, opts).map_err(|err| MicError::ExactFailed(err.to_string()))
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
    fn mic_default_options() {
        let opts = MicOptions::default();
        assert_eq!(opts.robust_mode, RobustMode::Filtered);
        assert!(!opts.use_bvh);
    }

    #[test]
    fn mic_exact_square_radius() {
        let poly = square(10.0);
        let opts = MicOptions {
            robust_mode: RobustMode::Filtered,
            use_bvh: false,
            use_parallel: false,
        };
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
        let opts = MicOptions {
            robust_mode: RobustMode::Filtered,
            use_bvh: false,
            use_parallel: false,
        };
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
        let opts = MicOptions {
            robust_mode: RobustMode::Filtered,
            use_bvh: false,
            use_parallel: false,
        };
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
}
