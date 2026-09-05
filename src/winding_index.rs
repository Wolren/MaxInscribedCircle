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
    pub fn new() -> Self {
        WindingIndex
    }

    pub fn from_polygon(_poly: &Polygon<f64>) -> Self {
        WindingIndex
    }

    pub fn contains(&self, _x: f64, _y: f64) -> bool {
        true
    }
}
