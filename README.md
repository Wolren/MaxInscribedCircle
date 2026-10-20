<div align="center">

# MaxInscribedCircle

Finds the largest circle that fits inside a polygon. Holes act as boundary, so
the circle must clear them too.

![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-27%20passing-brightgreen.svg)](https://github.com/Wolren/MaxInscribedCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-8.4x%20faster-blue.svg)](https://github.com/Wolren/MaxInscribedCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxInscribedCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

The maximum inscribed circle is the largest circle that fits inside a polygon.
Every boundary segment limits it, so an interior ring limits it as well.

## API

| Function | Use |
|----------|-----|
| `maximum_inscribed_circle` | Solve one polygon |
| `maximum_inscribed_circle_multipolygon` | Keep the best component |
| `maximum_inscribed_circle_slice` | Solve a batch with rayon |
| `maximum_inscribed_circle_with_workspace` | Reuse a built workspace |

`MicResult` reports the centre and the radius. It also reports the squared
radius, the indices of the segments the circle touches, and the candidate count.
`used_engine` and `approach` name the stage that answered. `component_index`
names the multipolygon component.

The batch entry point needs no feature, because rayon is a plain dependency.

Run `cargo run --example quickstart` for a complete program.

## How it solves

The crate runs two stages.

The first stage handles a single ring of four or five vertices. It solves the
triangle incenter, the axis-aligned rectangle, and the convex quad in closed
form. It then measures each result against the segment table. A result that
misses by more than 1e-9 relative falls through to the second stage.

The second stage is a grid solver. It runs branch and bound over the bounding
box and orders cells by the 1-Lipschitz bound on boundary distance. It refines
to `diag * 1e-7`. A compass hill climb then sharpens the incumbent. The search
pops at most 50,000 cells, and the climb finishes a long thin shape that the
budget cannot resolve.

The reported radius is always measured, never inferred. The crate measures the
distance from the centre to the nearest boundary segment at the end of the
solve, so a returned circle fits the polygon.

## Correctness

`tests/corpus_invariants.rs` checks two properties.

The first property is fit. The reported radius must equal the measured clearance
at the reported centre.

The second property is reach. No dense interior sample may clear more than the
returned circle.

The workflow in `.github/workflows/ci.yml` runs the suite on every push and pull
request. It runs four feature sets, the 1.88 MSRV check, and a clippy, fmt, and
rustdoc job. No check needs an external solver.

## Performance

The corpus holds 290 Polish cadastral parcels in
`benches/data/realworld_290.geojson`. Both solvers receive the same polygon and
the same tolerance. The timer covers the solve only, because both sides build
their indexes before the clock starts. Ratio is GEOS divided by native. Both
sides run serially, so the table measures the solver and not the threading.

| Measure | native | GEOS | Ratio |
|---------|-------:|-----:|------:|
| total, 290 parcels | **41.4 ms** | 347.5 ms | **8.39x** |
| per-case median | | | 6.3x |
| per-case p10 | | | 2.0x |
| per-case p90 | | | 14.2x |
| native faster on | 279/290 | | 96% |

The margin narrows as vertex count grows.

| Exterior vertices | n | native | GEOS | Ratio |
|------------------:|--:|-------:|-----:|------:|
| 0-7 | 42 | 1.99 ms | 30.6 ms | 15.3x |
| 8-16 | 111 | 16.7 ms | 226 ms | 13.5x |
| 17-32 | 75 | 8.95 ms | 74.6 ms | 8.34x |
| 33-64 | 28 | 4.01 ms | 13.5 ms | 3.36x |
| 129+ | 14 | 6.91 ms | 13.1 ms | 1.90x |

Pinned to one core the ratio holds at 7.9x. The batch entry point is separate.
`maximum_inscribed_circle_slice` returns bit-identical results in input order and
runs 6.1x faster on 12 cores.

Repeat runs of the same tree land between 8.39x and 8.78x, so the second digit is
noise. Regenerate before quoting a new figure.

Regenerate with
`cargo run --release --features bench-geos --example compare_apples`. The build
needs the GEOS C library. The GEOS dependency reads `GEOS_LIB_DIR` and
`GEOS_VERSION` when pkg-config cannot find an install.

## Integration with the geo ecosystem

Enable `geo-traits` to accept a foreign polygon directly. WKB views, GeoArrow
arrays, and GeoZero streams then need no conversion.

| Function | Accepts |
|----------|---------|
| `maximum_inscribed_circle_polygon_trait` | any `PolygonTrait` with `f64` coordinates |
| `maximum_inscribed_circle_multipolygon_trait` | any `MultiPolygonTrait` with `f64` coordinates |
| `maximum_inscribed_circle_geometry` | `geo_types::Geometry`, polygon and multipolygon only |

The dispatch function rejects every other variant with `MicError::InvalidInput`
and names the variant in the message.

## Python bindings

The `python` feature builds the extension module `max_inscribed_circle`. Build
it with maturin, or with `cargo build --features python,extension-module`.

| Function | Signature |
|----------|-----------|
| `solve_mic_py` | `(exterior, holes=None)` |
| `solve_mic_multipolygon_py` | `(polygons)` |

The result class reports `center_x`, `center_y`, `radius`, `radius_sq`,
`used_engine`, `candidate_count`, and `component_index`. Rings arrive as
`(x, y)` tuples. A rejected input raises `ValueError`, and so does a solver
panic.

## C API

`src/ffi.rs` is always compiled. It exposes `mic_solve`,
`mic_solve_multipolygon`, and `mic_options_default`. Flat `[x0, y0, ...]`
coordinate arrays go in, and `#[repr(C)]` structs come out.

Each function returns 0 on success and -1 on failure. A panic becomes -1, and a
null or short buffer is rejected instead of crashing. `MicOptionsC` carries no
tunable field yet. It exists so the ABI has a stable options pointer.

The crate builds as `rlib`, `cdylib`, and `staticlib`, so `cargo build
--release` produces the shared library without an extra feature.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `bench-geos` | Comparison only: links the GEOS C library so benchmarks can cross-check the native solver. Never compiled into the solve path. | no |
| `geo-traits` | Generic solves over any `geo_traits` polygon, plus `geo_types::Geometry` dispatch | no |
| `shewchuk` | Exact `orient2d` predicates from `robust` in the containment and orientation tests | no |
| `simd` | 4-lane batched segment distance through the stable `wide` crate | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_inscribed_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |

## License

Apache-2.0
