//! Cross-check the native solvers against GEOS on a handful of shapes.
//!
//! Build the GEOS C library first, then run:
//!
//! ```text
//! GEOS_LIB_DIR=<lib dir> GEOS_INCLUDE_DIR=<include dir> GEOS_VERSION=3.14.1 \
//!   cargo run --release --features bench-geos --example compare_geos
//! ```
//!
//! Both solvers see the same coordinates and the report shows the radius each
//! returns, the relative gap between them, and the time each takes per solve.
//! GEOS runs the same cell-refinement branch and bound as the native grid
//! stage, so a small gap shows two implementations of one algorithm agree, not
//! that an independent method confirmed the radius.
//!
//! Timings are indicative: one warm-up call each, then a small fixed number of
//! repeats on an otherwise idle machine.

use std::time::Instant;

use geo_types::{Coord, LineString, Polygon};
use max_inscribed_circle::mic::input::HostPolygon;
use max_inscribed_circle::mic::solver::bench_geos::{geos_tolerance_for, solve_mic_geos};
use max_inscribed_circle::mic::{
    maximum_inscribed_circle, MicApproach, MicOptions, MicResult, MicUsedEngine,
};

const NATIVE_REPEATS: u32 = 200;
const GEOS_REPEATS: u32 = 20;

fn ring(points: &[(f64, f64)]) -> LineString<f64> {
    LineString::from(
        points
            .iter()
            .map(|&(x, y)| Coord { x, y })
            .collect::<Vec<_>>(),
    )
}

fn square(size: f64) -> Polygon<f64> {
    Polygon::new(
        ring(&[
            (0.0, 0.0),
            (size, 0.0),
            (size, size),
            (0.0, size),
            (0.0, 0.0),
        ]),
        vec![],
    )
}

fn l_shape() -> Polygon<f64> {
    Polygon::new(
        ring(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (2.0, 2.0),
            (2.0, 6.0),
            (0.0, 6.0),
            (0.0, 0.0),
        ]),
        vec![],
    )
}

fn square_with_hole() -> Polygon<f64> {
    Polygon::new(
        ring(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]),
        vec![ring(&[
            (3.0, 3.0),
            (7.0, 3.0),
            (7.0, 7.0),
            (3.0, 7.0),
            (3.0, 3.0),
        ])],
    )
}

fn staircase() -> Polygon<f64> {
    Polygon::new(
        ring(&[
            (0.0, 0.0),
            (9.0, 0.0),
            (9.0, 3.0),
            (6.0, 3.0),
            (6.0, 6.0),
            (3.0, 6.0),
            (3.0, 9.0),
            (0.0, 9.0),
            (0.0, 0.0),
        ]),
        vec![],
    )
}

fn sliver() -> Polygon<f64> {
    Polygon::new(
        ring(&[(0.0, 0.0), (20.0, 0.0), (20.0, 0.25), (0.0, 0.0)]),
        vec![],
    )
}

fn engine_name(engine: MicUsedEngine) -> &'static str {
    match engine {
        MicUsedEngine::Grid => "grid",
        MicUsedEngine::Geos => "geos",
    }
}

fn approach_name(approach: MicApproach) -> &'static str {
    match approach {
        MicApproach::FastTriangle => "FastTriangle",
        MicApproach::FastRectangle => "FastRectangle",
        MicApproach::FastConvexQuad => "FastConvexQuad",
        MicApproach::FastConvexN => "FastConvexN",
        MicApproach::GridSolver => "GridSolver",
        MicApproach::Geos => "Geos",
    }
}

fn label(r: &MicResult) -> String {
    format!(
        "{}/{}",
        engine_name(r.used_engine),
        approach_name(r.approach)
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let options = MicOptions::default();

    let shapes: Vec<(&str, Polygon<f64>)> = vec![
        ("square 10", square(10.0)),
        ("L shape", l_shape()),
        ("square hole", square_with_hole()),
        ("staircase", staircase()),
        ("sliver", sliver()),
    ];

    println!(
        "{:<13} {:>12} {:>12} {:>11} {:>10} {:>9} {:>9}  {}",
        "shape",
        "native radius",
        "geos radius",
        "rel gap",
        "native us",
        "geos us",
        "tol",
        "native path"
    );

    for (name, poly) in &shapes {
        let host = HostPolygon::from_polygon(poly)?;
        let tolerance = geos_tolerance_for(&host);

        let native = maximum_inscribed_circle(poly, &options)?;
        let geos_result = solve_mic_geos(poly)?;

        // Warm both paths before timing, so neither pays a first-call cost.
        let _ = maximum_inscribed_circle(poly, &options)?;
        let _ = solve_mic_geos(poly)?;

        let native_start = Instant::now();
        for _ in 0..NATIVE_REPEATS {
            let _ = maximum_inscribed_circle(poly, &options)?;
        }
        let native_us = native_start.elapsed().as_secs_f64() * 1e6 / f64::from(NATIVE_REPEATS);

        let geos_start = Instant::now();
        for _ in 0..GEOS_REPEATS {
            let _ = solve_mic_geos(poly)?;
        }
        let geos_us = geos_start.elapsed().as_secs_f64() * 1e6 / f64::from(GEOS_REPEATS);

        let gap = if geos_result.radius > 0.0 {
            (native.radius - geos_result.radius) / geos_result.radius
        } else {
            0.0
        };

        println!(
            "{:<13} {:>12.6} {:>12.6} {:>11.2e} {:>10.1} {:>9.1} {:>9.1e}  {}",
            name,
            native.radius,
            geos_result.radius,
            gap,
            native_us,
            geos_us,
            tolerance,
            label(&native)
        );

        if geos_result.used_engine != MicUsedEngine::Geos {
            eprintln!("{}: geos path reported {:?}", name, geos_result.used_engine);
        }
    }

    Ok(())
}
