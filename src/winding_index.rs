//! Winding-number index stub carried over from IGE
//! `ige-core::solvers::common::winding_index`.
//!
//! Only wired behind the `shewchuk` feature in `mic::index::pip`.
//! The default ray-casting path does not touch this type.

use geo_types::Polygon;

/// Winding number index for point-in-polygon tests.
#[derive(Debug, Clone)]
pub struct WindingIndex;

impl WindingIndex {
    /// Create an index that reports every point as inside.
    pub fn new() -> Self {
        WindingIndex
    }

    /// Create an index for a polygon; the stub keeps no geometry.
    pub fn from_polygon(_poly: &Polygon<f64>) -> Self {
        WindingIndex
    }

    /// Return true for any query point; the stub runs no containment test.
    pub fn contains(&self, _x: f64, _y: f64) -> bool {
        true
    }
}
