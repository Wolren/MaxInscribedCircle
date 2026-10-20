# max-inscribed-circle

The largest circle that fits inside a polygon, holes included. Pure Rust with
no external solver in the pipeline.

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

## Install

```toml
[dependencies]
max-inscribed-circle = "0.1"
```

Rust 1.88 or newer. Default features are empty; the crate links no C library
and needs no system dependency.

## How it solves

Two stages, both native Rust. `MicResult::used_engine` and
`MicResult::approach` report which one answered.

1. **Closed-form fast paths** for single-ring inputs of four or five vertices:
   triangle incenter, axis-aligned rectangle, then a convex-quad Chebyshev
   solve. Every fast-path radius is re-measured against the segment table
   before it returns, and a result that disagrees by more than 1e-9 relative
   falls through to the grid instead of being trusted.
2. **Grid solver**: branch and bound over the bounding box, ordered by an upper
   bound from the 1-Lipschitz property of boundary distance, refining to a
   tolerance of `diag * 1e-7`. A compass hill climb on the clearance field
   polishes the incumbent, which is what keeps long thin shapes accurate when
   the cell budget would otherwise stop the bisection a percent short.

The returned radius is always the true distance from the returned centre to
the nearest boundary segment, measured at the end of the solve. A circle this
crate returns fits inside the polygon; it is never a bound read off a cell.

## Why you can trust the radius

Two independent properties are checked on every push, without needing any
external solver, in `tests/corpus_invariants.rs`:

- **fits**: the reported radius equals the measured clearance at the reported
  centre, to 1e-9 relative.
- **not short**: no dense sample of interior points clears more than the
  returned circle.

A third harness, `examples/parity_geos.rs`, cross-checks against GEOS
`MaximumInscribedCircle` over a corpus of 906 polygons (real WKT fixtures plus
synthetic shapes). Last run on this tree:

```
relative gap vs GEOS     max 1.502e-4   p90 2.233e-7   median 2.352e-8
radius above true clearance           native 0   GEOS 0
returned centre outside the polygon   native 0   GEOS 0
sampled floor above a returned radius native 0   GEOS 0
wall time over the corpus   native 2.99s   GEOS 23.2s
```

That harness needs the `bench-geos` feature and a GEOS install. It is
comparison only: nothing in the solve pipeline calls GEOS. GEOS runs the same
cell-refinement branch and bound as the grid stage, so a close match means two
implementations of one algorithm agree, not that an independent method
confirms the radius. The independent evidence is the two invariants above.

## API

| Function | Use |
|----------|-----|
| `maximum_inscribed_circle` | one polygon |
| `maximum_inscribed_circle_multipolygon` | best circle across components |
| `maximum_inscribed_circle_slice` | rayon-parallel batch, input order preserved |
| `maximum_inscribed_circle_with_workspace` | reuse a `MicWorkspace` across calls |

`MicResult` carries the centre, the radius and its square, the indices of the
supporting segments, the cell count the search popped, the solver and approach
that answered, and the component index for multipolygon solves.

## Features

| Feature | Effect |
|---------|--------|
| `geo-traits` | solve over any `geo_traits::PolygonTrait` or `MultiPolygonTrait`, plus `geo_types::Geometry` dispatch |
| `shewchuk` | exact `orient2d` predicates from `robust` in the containment test |
| `simd` | 4-lane batched segment distance over the stable `wide` crate |
| `hotpath` | `#[hotpath::measure]` instrumentation on the solver entry points |
| `python` | PyO3 extension module |
| `extension-module` | build the Python module as an importable extension |
| `bench-geos` | comparison only: links GEOS so benchmarks can cross-check the solvers |

`shewchuk` resolves containment queries that land within rounding distance of
an edge by exact arithmetic instead of by the sign of a cancellation error.
`simd` batches the 4-wide segment distance inner loop; it is a pure speed
feature and produces the same answers.

## Geo-traits interop

With `geo-traits` the solvers accept any foreign polygon without materialising
a `geo_types::Polygon` first, so WKB views, GeoArrow arrays and GeoZero streams
go straight in:

```rust
use max_inscribed_circle::geo_interop::maximum_inscribed_circle_polygon_trait;
use max_inscribed_circle::MicOptions;

let mic = maximum_inscribed_circle_polygon_trait(&foreign_polygon, &MicOptions::default())?;
```

`maximum_inscribed_circle_geometry` dispatches a `geo_types::Geometry`:
`Polygon` solves directly, `MultiPolygon` keeps the best component, and every
other variant is rejected with an error naming it rather than silently
mis-solved.

## Other interfaces

`mic_solve`, `mic_solve_multipolygon` and `mic_options_default` form a C ABI in
`src/ffi.rs`, and the crate builds as `rlib`, `cdylib` and `staticlib`. Panics
are caught at the boundary and become `-1`. The `python` feature exposes
`solve_mic_py` and `solve_mic_multipolygon_py`.

## Limits

- Input must be a valid simple polygon. Self-intersecting rings are not
  detected and are not repaired.
- Coordinates must be finite. A NaN or infinite coordinate is rejected with
  `MicError::InvalidInput` instead of returning a meaningless circle.
- A degenerate ring, such as zero area, is rejected the same way.
- The search is bounded: a shape whose optimum is a long ridge can spend its
  whole cell budget before the bisection closes, and the hill climb then
  finishes the job. Such a shape costs more than a compact one of the same
  vertex count. Accuracy stays within the tolerance above; on the GEOS corpus
  the worst case is 1.5e-4 relative.
- Accuracy is relative to the size of the polygon, not absolute. A polygon
  stored in degrees refines to a tolerance in degrees, so the same shape in
  metres and in degrees converges to the same circle.

## Examples

`cargo run --example quickstart` solves a square and a multipolygon.
`examples/probe_cost.rs` reports per-stage cost across a shape sweep.
`examples/parity_geos.rs` runs the GEOS comparison; set `MIC_WKT_DIRS` or pass
fixture directories as arguments.

## License

Apache-2.0
