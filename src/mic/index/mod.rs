//! Spatial indexes for the MIC solvers: point-in-polygon tests and
//! nearest-boundary distance queries.
pub mod bvh;
pub mod nearest_boundary;
pub mod pip;

pub use nearest_boundary::NearestBoundaryIndex;
pub use pip::PipIndex;
