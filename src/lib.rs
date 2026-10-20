//! Maximum inscribed circle (MIC) solvers for polygons.
//!
//! The maximum inscribed circle is the largest circle that fits entirely inside
//! a polygon. Holes are respected: the circle has to clear interior rings as
//! well as the exterior, so the solve runs against every boundary segment.
//!
//! # Quick start
//!
//! ```
//! use max_inscribed_circle::geo_types::{LineString, Polygon};
//! use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};
//!
//! let square = Polygon::new(
//!     LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
//!     vec![],
//! );
//!
//! let mic = maximum_inscribed_circle(&square, &MicOptions::default())?;
//! assert!((mic.center.x() - 2.0).abs() < 1e-9);
//! assert!((mic.radius - 2.0).abs() < 1e-9);
//! # Ok::<(), max_inscribed_circle::MicError>(())
//! ```
//!
//! # Solve pipeline
//!
//! Three stages, all native Rust. [`MicResult::used_engine`] reports which one
//! produced the answer:
//!
//! 1. **Analytical fast paths** for single-ring inputs of four or five
//!    vertices: triangle incenter, then rectangle, then a convex quad solve.
//!    The convex quad result is checked against true edge distances and
//!    discarded when it is more than 1% off.
//! 2. **Grid search** over the bounding box, with a tolerance derived from the
//!    diagonal (`max(diag * 1e-7, 1e-12)`). It returns early when it finds a
//!    candidate, which is the usual outcome for complex shapes.
//! 3. **Exact Voronoi solve** over candidate sites, reached when the grid stage
//!    finds nothing.
//!
//! A degenerate input, such as a zero-area ring, fails fast with
//! [`MicError::InvalidInput`] instead of returning a meaningless circle.
//!
//! # Entry points
//!
//! | Function | Use |
//! |----------|-----|
//! | [`maximum_inscribed_circle`] | one polygon |
//! | [`maximum_inscribed_circle_multipolygon`] | best circle across components |
//! | [`maximum_inscribed_circle_slice`] | rayon-parallel batch, input order preserved |
//! | [`maximum_inscribed_circle_with_workspace`] | reuse a [`MicWorkspace`] across calls |
//!
//! # GEOS comparison
//!
//! The `bench-geos` feature compiles `mic::solver::bench_geos`, which calls
//! the GEOS `MaximumInscribedCircle` construct so benchmarks can cross-check
//! the native solvers. It is not an engine: nothing in [`mic`] calls it, and
//! the default feature set does not compile it. GEOS runs the same
//! cell-refinement branch and bound as the grid stage, so agreement between
//! them means two implementations of one algorithm match, not that an
//! independent method confirms the radius.
//!
//! # Features
//!
//! | Feature | Effect |
//! |---------|--------|
//! | `bench-geos` | comparison only: links the GEOS C library so benchmarks can compare solvers |
//! | `geo-traits` | generic solves over any `geo_traits` polygon, plus `geo_types::Geometry` dispatch |
//! | `shewchuk` | exact `orient2d` predicates from the `robust` crate in the containment and orientation tests |
//! | `simd` | 4-lane batched segment distance via the stable `wide` crate |
//! | `hotpath` | `#[hotpath::measure]` instrumentation on the solver entry points |
//! | `python` | PyO3 extension module named `max_inscribed_circle` |
//! | `extension-module` | build the Python module as an importable extension |
//!
//! Default features are empty. The `bench-geos` feature links against the GEOS
//! C library, and its build script expects `GEOS_LIB_DIR` and `GEOS_VERSION`
//! when GEOS cannot be found through `pkg-config`.
//!
//! # C ABI
//!
//! [`ffi`] exposes `mic_solve`, `mic_solve_multipolygon` and
//! `mic_options_default`: flat coordinate arrays in, `#[repr(C)]` structs out,
//! `0` on success and `-1` on error, with panics caught at the boundary. The
//! crate builds as `rlib`, `cdylib` and `staticlib`.

#![warn(missing_docs)]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod ffi;
pub mod mic;
pub mod tuning;

#[cfg(feature = "python")]
#[cfg_attr(docsrs, doc(cfg(feature = "python")))]
pub mod python;

#[cfg(feature = "geo-traits")]
#[cfg_attr(docsrs, doc(cfg(feature = "geo-traits")))]
pub mod geo_interop;

/// Pinned `geo-types` re-export so downstream crates lock to the same
/// version this crate solves over.
pub use geo_types;

pub use mic::workspace::MicWorkspace;
pub use mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon,
    maximum_inscribed_circle_slice, maximum_inscribed_circle_with_workspace, MicApproach, MicError,
    MicOptions, MicResult, MicUsedEngine, RobustMode,
};

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Python extension module (feature `python`).
///
/// Build with maturin or with
/// `cargo build --features python,extension-module`, then import
/// `max_inscribed_circle`.
#[cfg(feature = "python")]
#[pymodule]
pub fn max_inscribed_circle(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    python::register(py, m)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
