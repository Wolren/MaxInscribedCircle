//! MIC solver backends: the closed-form fast paths and the grid solver.
//!
//! `bench_geos` sits alongside them behind the `bench-geos` feature, for
//! benchmark and comparison runs only; the solve pipeline never calls it.
pub mod fast;
pub mod grid;

#[cfg(feature = "bench-geos")]
pub mod bench_geos;
