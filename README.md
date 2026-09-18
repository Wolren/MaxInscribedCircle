# max-inscribed-circle

Largest circle that fits inside a polygon, holes included. Pure Rust, no
external solver.

```rust
use max_inscribed_circle::geo_types::{LineString, Polygon};
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

fn main() -> Result<(), max_inscribed_circle::MicError> {
    let square = Polygon::new(
        LineString::from(vec![(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0), (0.0, 0.0)]),
        vec![],
    );

    let mic = maximum_inscribed_circle(&square, &MicOptions::default())?;
    assert!((mic.center.x() - 2.0).abs() < 1e-9);
    assert!((mic.radius - 2.0).abs() < 1e-9);
    Ok(())
}
```

The solve runs analytical fast paths for triangles, rectangles and convex
quads, then a grid search over the bounding box, then an exact candidate
enumeration over Voronoi sites. `MicResult::used_engine` reports which stage
answered. The crate documentation describes the full pipeline.

## API

| Function | Use |
|----------|-----|
| `maximum_inscribed_circle` | one polygon |
| `maximum_inscribed_circle_multipolygon` | best circle across components |
| `maximum_inscribed_circle_slice` | rayon-parallel batch, input order preserved |
| `maximum_inscribed_circle_with_workspace` | reuse a `MicWorkspace` across calls |

`MicResult` carries the center, the radius and its square, the indices of the
supporting edges, the candidate count, the solver and approach that answered,
and the component index for multipolygon solves.

## Features

| Feature | Effect |
|---------|--------|
| `bench-geos` | comparison only: links the GEOS C library so benchmarks can compare solvers |
| `geo-traits` | generic solves over any `geo_traits` polygon, plus `geo_types::Geometry` dispatch |
| `shewchuk` | reserved: swaps the containment index for a stub that reports every query point as inside |
| `python` | PyO3 extension module |

Default features are empty, and the published crate never links a C library.
`bench-geos` is for benchmark and comparison runs: it compiles
`mic::solver::bench_geos`, which calls GEOS `MaximumInscribedCircle` so a
harness can cross-check the native grid stage. Nothing in the solve pipeline
calls it. GEOS runs the same cell-refinement branch and bound as that grid
stage, so a close match means two implementations of one algorithm agree.
The `shewchuk` feature is carried over from IGE and is not usable yet.

## Other interfaces

`mic_solve`, `mic_solve_multipolygon` and `mic_options_default` form a C ABI in
`src/ffi.rs`, and the crate builds as `rlib`, `cdylib` and `staticlib`. The
`python` feature exposes `solve_mic_py` and `solve_mic_multipolygon_py`.

## Examples

`examples/quickstart.rs` solves a square and a multipolygon, and prints the
solver and approach behind each answer.

## Provenance

The solver sources come from the IGE project, split out here as a standalone
crate.

## License

Apache-2.0
