//! MIC solver backends: the exact candidate solver and the grid solver.
//!
//! `bench_geos` sits alongside them behind the `bench-geos` feature, for
//! benchmark and comparison runs only; the solve pipeline never calls it.
pub mod exact;
pub mod grid;

#[cfg(feature = "bench-geos")]
pub mod bench_geos;
