//! Maximum inscribed circle (MIC) solvers extracted from IGE.
//!
//! Pipeline: analytical fast paths (triangle, rectangle, convex quad) first,
//! then the grid solver, then the Voronoi exact solver over a constrained
//! Delaunay triangulation, with an optional GEOS fallback behind the `geos`
//! feature.
//!
//! ## Entry points
//! - `maximum_inscribed_circle` - solve on a single polygon
//! - `maximum_inscribed_circle_multipolygon` - best circle across components
//! - `maximum_inscribed_circle_slice` - rayon-parallel batch solve
//! - `maximum_inscribed_circle_with_workspace` - reuse a `MicWorkspace`
//!
//! ## Features
//! | Feature | Effect |
//! |---------|--------|
//! | `geos` | GEOS fallback engine |
//! | `shewchuk` | Winding-number point-in-polygon index |
//! | `python` | PyO3 bindings module `max_inscribed_circle` |
//! | `extension-module` | Build the Python module as an importable extension |
//! | `geo-traits` | Generic solves over `geo_traits` polygons plus `Geometry` dispatch |
//!
//! ## Foreign-function interface
//!
//! The C ABI lives in [`ffi`] (`mic_solve`, `mic_solve_multipolygon`):
//! flat coordinate arrays in, `#[repr(C)]` result structs out, 0 on
//! success, -1 on error, panics caught at the boundary.

pub mod ffi;
pub mod mic;
pub mod tuning;
pub mod winding_index;

#[cfg(feature = "python")]
pub mod python;

#[cfg(feature = "geo-traits")]
pub mod geo_interop;

/// Pinned `geo-types` re-export so downstream crates lock to the same
/// version this crate solves over.
pub use geo_types;

pub use mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon,
    maximum_inscribed_circle_slice, maximum_inscribed_circle_with_workspace, MicApproach, MicEngine,
    MicError, MicOptions, MicResult, MicUsedEngine, RobustMode,
};
pub use mic::workspace::MicWorkspace;

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
