//! MIC solver backends: the exact candidate solver, the grid solver, and
//! the optional GEOS fallback.
pub mod exact;
pub mod grid;

#[cfg(feature = "geos")]
pub mod geos_fallback;
