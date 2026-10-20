//! Real-corpus benchmark: characterise the 290 parcels, then time the solver.
//!
//! Run: cargo run --release --example real_corpus
//!
//! The first table is the important one. Optimising before knowing the shape
//! distribution is how a change that helps one polygon class quietly hurts the
//! class that actually dominates the corpus.

use std::time::Instant;

use max_inscribed_circle::mic::grid_tolerance;
use max_inscribed_circle::mic::input::HostPolygon;
use max_inscribed_circle::mic::solver::grid::solve_grid;
use max_inscribed_circle::mic::workspace::MicWorkspace;
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

#[path = "../benches/corpus.rs"]
mod corpus;

fn pct(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let i = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[i.min(sorted.len() - 1)]
}

fn main() {
    let parcels = corpus::load_parcels();
    println!("parcels: {}", parcels.len());

    // ---- shape distribution ----
    let mut verts: Vec<f64> = parcels
        .iter()
        .map(|p| corpus::exterior_len(&p.polygon) as f64)
        .collect();
    verts.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let holes = parcels
        .iter()
        .filter(|p| !p.polygon.interiors().is_empty())
        .count();

    let mut spans: Vec<f64> = parcels
        .iter()
        .map(|p| {
            let xs: Vec<f64> = p.polygon.exterior().0.iter().map(|c| c.x).collect();
            let ys: Vec<f64> = p.polygon.exterior().0.iter().map(|c| c.y).collect();
            let w = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - xs.iter().cloned().fold(f64::INFINITY, f64::min);
            let h = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - ys.iter().cloned().fold(f64::INFINITY, f64::min);
            w.max(h)
        })
        .collect();
    spans.sort_by(|a, b| a.partial_cmp(b).unwrap());

    println!("\n=== shape distribution ===");
    println!(
        "exterior vertices   min {}  p50 {}  p90 {}  p99 {}  max {}",
        verts[0],
        pct(&verts, 0.50),
        pct(&verts, 0.90),
        pct(&verts, 0.99),
        verts[verts.len() - 1]
    );
    println!(
        "bounding span (deg) min {:.6}  p50 {:.6}  p90 {:.6}  max {:.6}",
        spans[0],
        pct(&spans, 0.50),
        pct(&spans, 0.90),
        spans[spans.len() - 1]
    );
    println!("parcels with holes: {holes}");

    // bucket by vertex count
    let buckets = [
        (0usize, 8usize),
        (8, 16),
        (16, 32),
        (32, 64),
        (64, 128),
        (128, 100_000),
    ];
    println!(
        "\n{:>10} {:>6} {:>12} {:>12}",
        "verts", "n", "median us", "max us"
    );
    for (lo, hi) in buckets {
        let mut us: Vec<f64> = Vec::new();
        for p in &parcels {
            let v = corpus::exterior_len(&p.polygon);
            if v < lo || v >= hi {
                continue;
            }
            let t = Instant::now();
            let _ = maximum_inscribed_circle(&p.polygon, &MicOptions::default());
            us.push(t.elapsed().as_secs_f64() * 1e6);
        }
        if us.is_empty() {
            continue;
        }
        us.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "{:>10} {:>6} {:>12.1} {:>12.1}",
            format!("{lo}-{hi}"),
            us.len(),
            pct(&us, 0.5),
            us[us.len() - 1]
        );
    }

    // ---- whole-corpus timing ----
    let t = Instant::now();
    let mut ok = 0usize;
    let mut err = 0usize;
    let mut cells_total = 0usize;
    for p in &parcels {
        match maximum_inscribed_circle(&p.polygon, &MicOptions::default()) {
            Ok(r) => {
                ok += 1;
                cells_total += r.candidate_count;
            }
            Err(_) => err += 1,
        }
    }
    let total = t.elapsed().as_secs_f64();
    println!("\n=== whole corpus ===");
    println!("solved {ok}  failed {err}   total {total:.3} s   cells popped {cells_total}");

    // ---- stage split, to see where time goes ----
    let mut t_fast = 0.0f64;
    let mut t_grid_only = 0.0f64;
    let mut t_build = 0.0f64;
    let mut stages = std::collections::BTreeMap::new();
    for p in &parcels {
        let t = Instant::now();
        if let Ok(r) = maximum_inscribed_circle(&p.polygon, &MicOptions::default()) {
            *stages.entry(format!("{:?}", r.approach)).or_insert(0usize) += 1;
        }
        t_fast += t.elapsed().as_secs_f64();

        let Ok(host) = HostPolygon::from_polygon(&p.polygon) else {
            continue;
        };
        let t = Instant::now();
        let Ok(ws) = MicWorkspace::new(host.clone()) else {
            continue;
        };
        t_build += t.elapsed().as_secs_f64();

        let tol = grid_tolerance(&host);
        let t = Instant::now();
        let _ = solve_grid(&host, &ws.pip_index, &ws.nb_index, tol);
        t_grid_only += t.elapsed().as_secs_f64();
    }
    println!("\nfull pipeline {t_fast:.3} s   index build {t_build:.3} s   grid alone {t_grid_only:.3} s");
    println!("stages: {stages:?}");
}
