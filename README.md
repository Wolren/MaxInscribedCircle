<div align="center">

# MaxInscribedCircle

Finds the largest circle that fits inside a polygon.

[![crate](https://img.shields.io/crates/v/max-inscribed-circle.svg)](https://crates.io/crates/max-inscribed-circle)
[![docs](https://docs.rs/max-inscribed-circle/badge.svg)](https://docs.rs/max-inscribed-circle)
![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-27%20passing-brightgreen.svg)](https://github.com/Wolren/MaxInscribedCircle/actions/workflows/ci.yml)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-8.4x%20faster-blue.svg)](https://github.com/Wolren/MaxInscribedCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxInscribedCircle?tab=License-1-ov-file)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

The maximum inscribed circle is the largest circle that fits inside a polygon.
Holes act as boundary, so the circle must clear interior rings as well as the
exterior. This is the dual of the largest empty circle, which sits between
obstacles instead of inside them. The reported radius is measured at the returned
centre, so a returned circle fits the polygon. A zero-area ring or a non-finite
coordinate fails fast instead of returning a meaningless circle. The full API
reference and the algorithm detail are in the crate documentation:
`cargo doc --open`, no feature needed.

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
their indexes before the clock starts. Ratio is GEOS divided by native, so above 1
means this crate wins. Both sides run serially, so the table measures the solver
and not the threading.

| Measure | MaxInscribedCircle | GEOS | Ratio |
|---------|-------------------:|-----:|------:|
| total, 290 parcels | **41.4 ms** | 347.5 ms | **8.39x** |
| cases native solves faster | 279 of 290 | 11 of 290 | 96% |

The ratio varies by case. Across the 290 parcels:

| Per-case ratio | Value |
|----------------|------:|
| p10 | 2.0x |
| median | 6.3x |
| p90 | 14.2x |

The margin narrows as vertex count grows.

| Exterior vertices | n | MaxInscribedCircle | GEOS | Ratio |
|------------------:|--:|-------------------:|-----:|------:|
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

The `python` feature builds the extension module `max_inscribed_circle`. Build it
with maturin, or with `cargo build --features python,extension-module`.

| Function | Signature |
|----------|-----------|
| `solve_mic_py` | `(exterior, holes=None)` |
| `solve_mic_multipolygon_py` | `(polygons)` |

The result class reports `center_x`, `center_y`, `radius`, `radius_sq`,
`used_engine`, `candidate_count`, and `component_index`. Rings arrive as
`(x, y)` tuples. A rejected input raises `ValueError`, and so does a solver panic.

## C API

`src/ffi.rs` is always compiled, so the C ABI needs no feature. Build with
`cargo build --release`.

| Symbol | Signature |
|--------|-----------|
| `mic_solve` | `(coords, coords_len, holes, holes_len, options, out) -> c_int` |
| `mic_solve_multipolygon` | `(polys, polys_len, options, out) -> c_int` |
| `mic_options_default` | `() -> MicOptionsC` |

`coords` is flat `[x0, y0, ...]` exterior ring data. Holes arrive as an array of
`HoleInput`, each holding its own flat coordinate pointer. The multipolygon call
takes an array of `MicPolyInput`, each holding one exterior and its holes.
`MicResultC`, `MicOptionsC`, `HoleInput`, and `MicPolyInput` are `#[repr(C)]`.

Each function returns 0 on success and -1 on failure. A panic becomes -1, and a
null or short buffer is rejected instead of crashing.

## Features

| Feature | Description | Default |
|---------|-------------|---------|
| `geo-traits` | Generic solves over any `geo_traits` polygon, plus `geo_types::Geometry` dispatch | no |
| `simd` | 4-lane batched segment distance through the stable `wide` crate | no |
| `shewchuk` | Exact `orient2d` predicates from `robust` in the containment and orientation tests | no |
| `hotpath` | `hotpath::measure` instrumentation on the solver entry points | no |
| `python` | PyO3 extension module named `max_inscribed_circle` | no |
| `extension-module` | Build the Python module as an importable extension | no |
| `bench-geos` | Comparison only: links the GEOS C library so benchmarks can cross-check the native solver | no |

`rayon` is an unconditional dependency. The batch entry point
`maximum_inscribed_circle_slice` and the C ABI need no feature.

## License

Apache-2.0
