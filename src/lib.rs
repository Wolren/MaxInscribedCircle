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

pub mod mic;
pub mod tuning;
pub mod winding_index;

pub use mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon,
    maximum_inscribed_circle_slice, maximum_inscribed_circle_with_workspace, MicApproach, MicEngine,
    MicError, MicOptions, MicResult, MicUsedEngine, RobustMode,
};
pub use mic::workspace::MicWorkspace;
