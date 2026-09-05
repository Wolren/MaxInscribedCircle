# max-inscribed-circle

Maximum inscribed circle solvers extracted from [IGE](https://github.com/GIS-Inscribed-Geometry/IGE) `ige-core/src/solvers/mic`.

Pipeline: fast paths (triangle, rectangle, convex quad), grid solver, Voronoi exact solver over a constrained Delaunay triangulation (`spade`), optional GEOS fallback behind the `geos` feature.

## API

- `maximum_inscribed_circle` - single polygon
- `maximum_inscribed_circle_multipolygon` - best circle across components
- `maximum_inscribed_circle_slice` - rayon-parallel batch solve
- `maximum_inscribed_circle_with_workspace` - reusable workspace

## C ABI

`src/ffi.rs` exposes `mic_solve`, `mic_solve_multipolygon`, and
`mic_options_default`. Flat coordinate arrays in, `#[repr(C)]` structs
out, 0 on success, -1 on error, panics caught at the boundary. The
crate builds `cdylib` + `staticlib` + `rlib`. `FallbackOnly` without
the `geos` feature fails cleanly with -1.

## Python

Feature `python` provides the `max_inscribed_circle` module with
`solve_mic_py` and `solve_mic_multipolygon_py`. Build an importable
module with `cargo build --features python,extension-module` (or
maturin).

## Geo conventions

The Rust API solves over `geo_types::Polygon` (re-exported as
`max_inscribed_circle::geo_types` for version lock-in). Feature
`geo-traits` adds generic solves over any foreign polygon implementing
`geo_traits::PolygonTrait` plus `Geometry` dispatch
(`maximum_inscribed_circle_geometry`): Polygon solves directly,
MultiPolygon runs the native multipolygon solve, anything else is
rejected with a named error.

## Provenance

Source files under `src/mic/` are byte copies of the IGE MIC solver at extraction time. `src/tuning.rs` carries only the `MIC_*` constants. `src/winding_index.rs` carries the stub previously at `solvers::common::winding_index`, wired only behind the `shewchuk` feature.
