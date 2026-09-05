# ige-mic

Maximum inscribed circle solvers extracted from [IGE](https://github.com/GIS-Inscribed-Geometry/IGE) `ige-core/src/solvers/mic`.

Pipeline: fast paths (triangle, rectangle, convex quad), grid solver, Voronoi exact solver over a constrained Delaunay triangulation (`spade`), optional GEOS fallback behind the `geos` feature.

## API

- `maximum_inscribed_circle` - single polygon
- `maximum_inscribed_circle_multipolygon` - best circle across components
- `maximum_inscribed_circle_slice` - rayon-parallel batch solve
- `maximum_inscribed_circle_with_workspace` - reusable workspace

## Provenance

Source files under `src/mic/` are byte copies of the IGE MIC solver at extraction time. `src/tuning.rs` carries only the `MIC_*` constants. `src/winding_index.rs` carries the stub previously at `solvers::common::winding_index`, wired only behind the `shewchuk` feature.
