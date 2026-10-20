<div align="center">

# MaxInscribedCircle

Finds the largest circle that fits inside a polygon. Holes act as boundary, so
the circle must clear them too.

![MSRV](https://img.shields.io/badge/rustc-1.88+-ab6000.svg)
[![Tests](https://img.shields.io/badge/tests-27%20passing-brightgreen.svg)](https://github.com/Wolren/MaxInscribedCircle)
[![GEOS](https://img.shields.io/badge/vs%20GEOS-8.56x%20faster-blue.svg)](https://github.com/Wolren/MaxInscribedCircle#performance)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/Wolren/MaxInscribedCircle?tab=License-1-ov-file)
[![Status](https://img.shields.io/badge/status-experimental-orange.svg)]()
[![Ko-fi](https://img.shields.io/badge/Ko--fi-Support%20Wolren-FF5E5B?logo=ko-fi&logoColor=white)](https://ko-fi.com/wolren)

</div>

Pure Rust. The crate calls no external solver.

## Install

Add the dependency:

```toml
[dependencies]
max-inscribed-circle = "0.1"
```

Rust 1.88 or newer. Default features are empty. The crate links no C library.

## API

| Function | Use |
|----------|-----|
| `maximum_inscribed_circle` | Solve one polygon |
| `maximum_inscribed_circle_multipolygon` | Keep the best component |
| `maximum_inscribed_circle_slice` | Solve a batch with rayon |
| `maximum_inscribed_circle_with_workspace` | Reuse a built workspace |

`MicResult` reports the centre, the radius, the supporting segments, and the
stage that answered. A caller reads `used_engine` to learn which stage ran.

Run `cargo run --example quickstart` for a complete program.

## How it solves

The crate runs two stages.

The first stage handles single rings of four or five vertices. It solves the
triangle incenter, the axis-aligned rectangle, and the convex quad in closed
form. It then measures each result against the segment table. A result that
misses by more than 1e-9 relative falls through to the second stage.

The second stage is a grid solver. It runs branch and bound over the bounding
box and orders cells by the 1-Lipschitz bound on boundary distance. It refines
to `diag * 1e-7`. A compass hill climb then sharpens the incumbent. The climb is
what keeps a long thin shape accurate when the cell budget stops the bisection
early.

The reported radius is always measured, never inferred. The crate measures the
distance from the centre to the nearest boundary segment at the end of the
solve. It never reports a bound read from a cell.

## Correctness

`tests/corpus_invariants.rs` checks two properties on every push.

The first property is fit. The reported radius must equal the measured clearance
at the reported centre.

The second property is reach. No dense interior sample may clear more than the
returned circle.

Both checks run without an external solver.

## Performance

The corpus holds 290 Polish cadastral parcels in
`benches/data/realworld_290.geojson`. Both solvers receive the same polygon and
the same tolerance. The timer covers the solve only, because both sides build
their indexes before the clock starts. Ratio is GEOS divided by native.

| Measure | native | GEOS | Ratio |
|---------|-------:|-----:|------:|
| total, 290 parcels | **42.3 ms** | 361.8 ms | **8.56x** |
| per-case median | | | 6.6x |
| per-case p10 | | | 2.1x |
| native faster on | 280/290 | | 97% |

The margin narrows as vertex count grows.

| Exterior vertices | n | native | GEOS | Ratio |
|------------------:|--:|-------:|-----:|------:|
| 0-7 | 42 | 1.96 ms | 30.6 ms | 15.6x |
| 8-16 | 111 | 16.4 ms | 215 ms | 13.2x |
| 17-32 | 75 | 8.84 ms | 73.9 ms | 8.35x |
| 33-64 | 28 | 3.94 ms | 12.9 ms | 3.27x |
| 129+ | 14 | 6.94 ms | 12.8 ms | 1.84x |

Both sides run serially, so the table measures the solver and not the threading.
Pinned to one core the ratio holds at 7.9x. The batch entry point is separate.
`maximum_inscribed_circle_slice` returns bit-identical results in input order and
runs 6.1x faster on 12 cores.

Regenerate with `cargo run --release --features bench-geos --example compare_apples`.

## Features

| Feature | Effect |
|---------|--------|
| `geo-traits` | Solve over any `PolygonTrait` or `MultiPolygonTrait` |
| `shewchuk` | Use exact `orient2d` predicates from `robust` |
| `simd` | Batch segment distance over the `wide` crate |
| `hotpath` | Instrument the entry points with `hotpath` |
| `python` | Build a PyO3 extension module |
| `bench-geos` | Link GEOS for benchmark comparison only |

Enable `geo-traits` to accept a foreign polygon directly. WKB views, GeoArrow
arrays, and GeoZero streams then need no conversion. The C ABI lives in
`src/ffi.rs`, and it catches panic at the boundary.

## Limits

Input must be a valid simple polygon. The crate does not detect a self
intersection.

The crate rejects a NaN or infinite coordinate with `MicError::InvalidInput`. A
zero-area ring fails the same way.

A long ridge can exhaust the cell budget. The hill climb then finishes the work.
Accuracy stays within the tolerance above. The worst case on the GEOS corpus is
1.5e-4 relative.

Accuracy is relative to polygon size. A polygon stored in degrees refines to a
tolerance in degrees. The same shape in metres and in degrees agrees.

## License

Apache-2.0
