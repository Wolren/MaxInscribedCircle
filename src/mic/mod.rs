//! Maximum Inscribed Circle (MIC) solvers for polygonal inputs.

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

/// Engine selection for MIC solving.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MicEngine {
    /// Use the native Rust polygon-specialized solver only.
    ExactOnly,
    /// Use GEOS fallback only.
    FallbackOnly,
    /// Try exact first, then GEOS fallback if exact fails.
    ExactThenGeos,
}

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
    /// Engine selection: exact only, GEOS fallback only, or exact then GEOS.
    pub engine: MicEngine,
    /// Robustness mode: filtered adds candidate filtering and a tighter GEOS tolerance.
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
            engine: MicEngine::ExactThenGeos,
            robust_mode: RobustMode::Filtered,
            use_bvh: false,
            use_parallel: false,
        }
    }
}

/// Engine that produced the final result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MicUsedEngine {
    /// The native exact solver produced the circle.
    Exact,
    /// The grid solver produced the circle.
    Grid,
    /// The GEOS fallback produced the circle.
    GeosFallback,
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
    /// GEOS maximum inscribed circle fallback.
    GeosFallback,
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
    /// Engine that produced this result.
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
    /// The GEOS fallback was requested but the `geos` feature is disabled.
    #[error("GEOS fallback feature is not enabled")]
    GeosFeatureDisabled,
    /// The GEOS fallback failed.
    #[error("GEOS fallback failed: {0}")]
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
    if opts.engine != MicEngine::FallbackOnly && host.ring_count() == 1 {
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

    // Phase 1: Use grid solver for complex shapes - matches GEOS exactly and is much faster
    // Skip for FallbackOnly - let GEOS solve independently for comparison
    let use_grid = opts.engine != MicEngine::FallbackOnly;

    if use_grid {
        // Use grid solver - now uses original polygon (no normalization)
        let workspace = MicWorkspace::new(host.clone(), opts.use_bvh)?;
        let bounds = host.bounds().unwrap_or((0.0, 0.0, 1.0, 1.0));
        let diag = (bounds.2 - bounds.0).hypot(bounds.3 - bounds.1).max(1.0);
        let tolerance = (diag * 1e-7).max(1e-12);

        if let Some(grid_result) =
            solve_grid(host, &workspace.pip_index, &workspace.nb_index, tolerance)
        {
            return Ok(grid_result);
        }
    }

    match opts.engine {
        MicEngine::ExactOnly => run_exact(host, opts),
        MicEngine::FallbackOnly => run_geos(host, None, opts),
        MicEngine::ExactThenGeos => {
            let mut workspace = match MicWorkspace::new(host.clone(), opts.use_bvh) {
                Ok(w) => w,
                Err(e) => {
                    #[cfg(feature = "geos")]
                    {
                        return run_geos(host, None, opts).map_err(|fe| {
                            MicError::GeosFailed(format!(
                                "workspace failed ({e}); fallback failed ({fe})"
                            ))
                        });
                    }
                    #[cfg(not(feature = "geos"))]
                    {
                        return Err(e);
                    }
                }
            };
            match solve_exact(&mut workspace, opts) {
                Ok(result) => Ok(result),
                Err(e) => {
                    #[cfg(feature = "geos")]
                    {
                        let seg_index = workspace.nb_index.segments().clone();
                        run_geos(host, Some(&seg_index), opts).map_err(|fallback_err| {
                            MicError::GeosFailed(format!(
                                "exact failed ({e}); fallback failed ({fallback_err})"
                            ))
                        })
                    }
                    #[cfg(not(feature = "geos"))]
                    {
                        Err(MicError::ExactFailed(e.to_string()))
                    }
                }
            }
        }
    }
}

fn run_exact(host: &HostPolygon, opts: &MicOptions) -> Result<MicResult, MicError> {
    let mut workspace = MicWorkspace::new(host.clone(), opts.use_bvh)?;
    solve_exact(&mut workspace, opts).map_err(|err| MicError::ExactFailed(err.to_string()))
}

fn run_geos(
    host: &HostPolygon,
    existing_seg_index: Option<&input::SegmentIndex>,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    #[cfg(feature = "geos")]
    {
        self::solver::geos_fallback::solve_with_geos(host, opts, existing_seg_index)
    }
    #[cfg(not(feature = "geos"))]
    {
        let _ = (host, existing_seg_index, opts);
        Err(MicError::GeosFeatureDisabled)
    }
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
        assert_eq!(opts.engine, MicEngine::ExactThenGeos);
        assert_eq!(opts.robust_mode, RobustMode::Filtered);
        assert!(!opts.use_bvh);
    }

    #[test]
    fn mic_exact_square_radius() {
        let poly = square(10.0);
        let opts = MicOptions {
            engine: MicEngine::ExactOnly,
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
            engine: MicEngine::ExactOnly,
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
            engine: MicEngine::ExactOnly,
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
