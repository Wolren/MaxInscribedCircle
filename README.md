# max-inscribed-circle

Largest circle that fits inside a polygon, holes included. Pure Rust, with an
optional GEOS fallback.

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
quads, then a grid search, then an exact Voronoi solve over a constrained
Delaunay triangulation. `MicResult::used_engine` reports which stage answered.
The crate documentation describes the full pipeline.

## API

| Function | Use |
|----------|-----|
| `maximum_inscribed_circle` | one polygon |
| `maximum_inscribed_circle_multipolygon` | best circle across components |
| `maximum_inscribed_circle_slice` | rayon-parallel batch, input order preserved |
| `maximum_inscribed_circle_with_workspace` | reuse a `MicWorkspace` across calls |

`MicResult` carries the center, the radius and its square, the indices of the
supporting edges, the candidate count, the engine and approach that answered,
and the component index for multipolygon solves.

## Features

| Feature | Effect |
|---------|--------|
| `geo-traits` | generic solves over any `geo_traits` polygon, plus `geo_types::Geometry` dispatch |
| `geos` | GEOS fallback engine |
| `shewchuk` | reserved: swaps the containment index for a stub that reports every query point as inside |
| `python` | PyO3 extension module |

Default features are empty. `geos` links against the GEOS library at build time.
The `shewchuk` feature is carried over from IGE and is not usable yet.

## Other interfaces

`mic_solve`, `mic_solve_multipolygon` and `mic_options_default` form a C ABI in
`src/ffi.rs`, and the crate builds as `rlib`, `cdylib` and `staticlib`. The
`python` feature exposes `solve_mic_py` and `solve_mic_multipolygon_py`.

## Examples

`examples/quickstart.rs` solves a square and a multipolygon, and prints the
engine and approach behind each answer.

## Provenance

The solver sources come from the IGE project, split out here as a standalone
crate.

## License

Apache-2.0
