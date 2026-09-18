//! Where the time goes on one polygon, and what a batch run looks like when
//! the caller uses the rayon entry point instead of a sequential loop.
//!
//! Run: cargo run --release --example probe_cost

use std::time::Instant;

use geo_types::{Coord, LineString, Polygon};

use max_inscribed_circle::mic::grid_tolerance;
use max_inscribed_circle::mic::input::HostPolygon;
use max_inscribed_circle::mic::solver::exact::solve_exact;
use max_inscribed_circle::mic::solver::grid::solve_grid;
use max_inscribed_circle::mic::workspace::MicWorkspace;
use max_inscribed_circle::mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_slice, MicOptions,
};

fn ring(pts: &[(f64, f64)]) -> LineString<f64> {
    let mut v: Vec<Coord<f64>> = pts.iter().map(|(x, y)| Coord { x: *x, y: *y }).collect();
    if v.first() != v.last() {
        v.push(v[0]);
    }
    LineString::from(v)
}

fn poly(pts: &[(f64, f64)]) -> Polygon<f64> {
    Polygon::new(ring(pts), vec![])
}

/// Rectilinear staircase: `steps` alternating runs of length `run`.
fn staircase(steps: usize, run: f64) -> Polygon<f64> {
    let mut pts = vec![(0.0, 0.0)];
    let (mut x, mut y) = (0.0, 0.0);
    for s in 0..steps {
        if s % 2 == 0 {
            y += run;
            pts.push((x, y));
            x += run;
            pts.push((x, y));
        } else {
            x += run;
            pts.push((x, y));
            y += run;
            pts.push((x, y));
        }
    }
    pts.push((x, 0.0));
    pts.push((0.0, 0.0));
    poly(&pts)
}

fn rect(w: f64, h: f64) -> Polygon<f64> {
    poly(&[(0.0, 0.0), (w, 0.0), (w, h), (0.0, h)])
}

/// A long corridor with extra vertices on one edge. The extra vertices keep it
/// off the fast paths (which want 4 or 5 ring coordinates), and its optimum is
/// a ridge: every point along the mid line has the same clearance.
fn notched_corridor(len: f64, width: f64) -> Polygon<f64> {
    poly(&[
        (0.0, 0.0),
        (len, 0.0),
        (len, width),
        (len * 0.75, width),
        (len * 0.75, width * 0.94),
        (len * 0.5, width * 0.94),
        (len * 0.5, width),
        (0.0, width),
    ])
}

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}

fn main() {
    let opts = MicOptions::default();
    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);

    let shapes: Vec<(String, Polygon<f64>)> = vec![
        ("square 100".to_string(), rect(100.0, 100.0)),
        ("corridor 100x2".to_string(), rect(100.0, 2.0)),
        ("corridor 100x0.2".to_string(), rect(100.0, 0.2)),
        ("stair 3 runs x1".to_string(), staircase(3, 1.0)),
        ("stair 8 runs x1".to_string(), staircase(8, 1.0)),
        ("stair 8 runs x2".to_string(), staircase(8, 2.0)),
        ("stair 20 runs x2".to_string(), staircase(20, 2.0)),
        (
            "corridor+notch 100".to_string(),
            notched_corridor(100.0, 2.0),
        ),
        (
            "corridor+notch 1000".to_string(),
            notched_corridor(1000.0, 2.0),
        ),
        (
            "corridor+notch 100x0.2".to_string(),
            notched_corridor(100.0, 0.2),
        ),
    ];

    println!("cores {threads}   tolerance = 1e-7 of the bounding box diagonal\n");
    println!(
        "{:<18} {:>11} {:>9} {:>9} {:>9} {:>10} {:>10} {:>11} {:>9}",
        "shape",
        "radius",
        "ws ms",
        "grid ms",
        "exact ms",
        "full ms",
        "cells grid",
        "cells full",
        "tol"
    );

    for (name, p) in &shapes {
        let host = HostPolygon::from_polygon(p).expect("host");
        let tol = grid_tolerance(&host);

        let t = Instant::now();
        let mut ws = MicWorkspace::new(host.clone(), false).expect("workspace");
        let ws_ms = ms(t);

        let t = Instant::now();
        let grid = solve_grid(&host, &ws.pip_index, &ws.nb_index, tol);
        let grid_ms = ms(t);

        let t = Instant::now();
        let exact = solve_exact(&mut ws, &opts);
        let exact_ms = ms(t);

        let t = Instant::now();
        let full = maximum_inscribed_circle(p, &opts).expect("full");
        let full_ms = ms(t);

        let grid_radius = grid.as_ref().map(|r| r.radius).unwrap_or(f64::NAN);
        let exact_radius = exact.as_ref().map(|r| r.radius).unwrap_or(f64::NAN);
        let cells_grid = grid.as_ref().map(|r| r.candidate_count).unwrap_or(0);

        println!(
            "{:<18} {:>11.7} {:>9.2} {:>9.2} {:>9.2} {:>10.2} {:>10} {:>11} {:>9.1e}",
            name,
            full.radius,
            ws_ms,
            grid_ms,
            exact_ms,
            full_ms,
            cells_grid,
            full.candidate_count,
            tol
        );
        println!(
            "{:<18} {:>11} {:>11.7} {:>11.7}  grid radius | exact radius | full radius {:?}/{:?}",
            "", "", grid_radius, exact_radius, full.used_engine, full.approach
        );
    }

    // Batch: the same work, once as a sequential loop, once through the rayon
    // entry point the crate exposes for it.
    let corpus: Vec<Polygon<f64>> = (0..96)
        .map(|i| match i % 4 {
            0 => staircase(8, 1.0),
            1 => staircase(20, 2.0),
            2 => rect(100.0, 2.0),
            _ => staircase(3, 1.0),
        })
        .collect();

    let t = Instant::now();
    let mut seq_radius = 0.0;
    for p in &corpus {
        seq_radius += maximum_inscribed_circle(p, &opts).expect("seq").radius;
    }
    let seq_ms = ms(t);

    let t = Instant::now();
    let batch = maximum_inscribed_circle_slice(&corpus, &opts);
    let par_ms = ms(t);

    let batch_radius: f64 = batch
        .iter()
        .map(|r| r.as_ref().map(|r| r.radius).unwrap_or(0.0))
        .sum();
    let same = (batch_radius - seq_radius).abs() <= 1e-12 * seq_radius.abs().max(1.0);

    println!(
        "\n=== batch of {} staircase and corridor polygons ===",
        corpus.len()
    );
    println!("sequential loop   {seq_ms:8.1} ms");
    println!(
        "rayon batch       {par_ms:8.1} ms   speedup {:.2}x on {threads} cores",
        seq_ms / par_ms
    );
    println!("identical results {same}");

    // How each stage scales with the vertex count. Regular n-gons have an
    // isolated optimum at the centre, so this measures the stages themselves
    // rather than the shape of the optimum.
    let regular = |n: usize, r: f64| -> Polygon<f64> {
        let pts: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                (r * a.cos(), r * a.sin())
            })
            .collect();
        poly(&pts)
    };

    println!("\n=== vertex count sweep, regular n-gons of radius 10 ===");
    println!(
        "{:>6} {:>12} {:>11} {:>10} {:>11} {:>10} {:>10} {:>11} {:>9}",
        "verts",
        "true radius",
        "ws ms",
        "exact ms",
        "exact rel",
        "grid ms",
        "grid cells",
        "full ms",
        "engine"
    );
    for n in [4usize, 8, 16, 32, 64, 128, 256, 512, 1024] {
        let p = regular(n, 10.0);
        let host = HostPolygon::from_polygon(&p).expect("host");
        let tol = grid_tolerance(&host);

        let t = Instant::now();
        let mut ws = MicWorkspace::new(host.clone(), false).expect("workspace");
        let ws_ms = ms(t);

        let t = Instant::now();
        let ex = solve_exact(&mut ws, &opts).ok();
        let exact_ms = ms(t);

        let t = Instant::now();
        let gr = solve_grid(&host, &ws.pip_index, &ws.nb_index, tol);
        let grid_ms = ms(t);

        let t = Instant::now();
        let full = maximum_inscribed_circle(&p, &opts).expect("full");
        let full_ms = ms(t);

        let truth = 10.0 * (std::f64::consts::PI / n as f64).cos();
        let exact_rel = ex
            .as_ref()
            .map(|r| (r.radius - truth) / truth)
            .unwrap_or(f64::NAN);
        println!(
            "{:>6} {:>12.9} {:>11.2} {:>10.2} {:>11.2e} {:>10.2} {:>10} {:>11.2} {:>9?}",
            n,
            truth,
            ws_ms,
            exact_ms,
            exact_rel,
            grid_ms,
            gr.as_ref().map(|r| r.candidate_count).unwrap_or(0),
            full_ms,
            full.used_engine
        );
    }
}
