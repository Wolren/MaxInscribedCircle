//! Apples-to-apples MIC benchmark against GEOS.
//!
//! Both sides get the same polygon and tolerance, and only the solve is timed.
//! The native side runs the real public entry point, so fast paths are included.
//! GEOS gets a geometry built before the clock starts. Neither side includes
//! coordinate conversion or output parsing.
//!
//! Run: cargo run --release --features bench-geos --example compare_apples

use std::time::Instant;

use geos::{Geom, Geometry};
use max_inscribed_circle::geo_types::Polygon;
use max_inscribed_circle::mic::grid_tolerance;
use max_inscribed_circle::mic::input::HostPolygon;
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

#[path = "../benches/corpus.rs"]
mod corpus;

fn to_wkt(poly: &Polygon<f64>) -> String {
    let mut s = String::from("POLYGON (");
    let mut first = true;
    for ring in std::iter::once(poly.exterior()).chain(poly.interiors()) {
        if !first {
            s.push_str(", ");
        }
        first = false;
        s.push('(');
        for (i, c) in ring.0.iter().enumerate() {
            if i > 0 {
                s.push_str(", ");
            }
            s.push_str(&format!("{} {}", c.x, c.y));
        }
        s.push(')');
    }
    s.push(')');
    s
}

struct Row {
    name: String,
    verts: usize,
    native_us: f64,
    geos_us: f64,
    r_native: f64,
    r_geos: f64,
}

fn pct(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let i = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[i.min(sorted.len() - 1)]
}

fn main() {
    let parcels = corpus::load_parcels();
    let opts = MicOptions::default();

    struct Case {
        name: String,
        verts: usize,
        poly: Polygon<f64>,
        tol: f64,
        geom: Geometry,
    }
    let mut cases: Vec<Case> = Vec::new();
    for p in &parcels {
        let Ok(host) = HostPolygon::from_polygon(&p.polygon) else {
            continue;
        };
        let tol = grid_tolerance(&host);
        let Ok(geom) = Geometry::new_from_wkt(&to_wkt(&p.polygon)) else {
            continue;
        };
        cases.push(Case {
            name: p.name.clone(),
            verts: corpus::exterior_len(&p.polygon),
            poly: p.polygon.clone(),
            tol,
            geom,
        });
    }
    println!("cases: {}\n", cases.len());

    for c in cases.iter().take(30) {
        let _ = maximum_inscribed_circle(&c.poly, &opts);
        let _ = c.geom.maximum_inscribed_circle(c.tol);
    }

    let mut rows: Vec<Row> = Vec::with_capacity(cases.len());
    for c in &cases {
        let mut native = f64::INFINITY;
        let mut geos = f64::INFINITY;
        let (mut rn, mut rg) = (0.0, 0.0);
        for _ in 0..5 {
            let t = Instant::now();
            let out = maximum_inscribed_circle(&c.poly, &opts);
            let us = t.elapsed().as_secs_f64() * 1e6;
            if us < native {
                native = us;
            }
            if let Ok(r) = out {
                rn = r.radius;
            }

            let t = Instant::now();
            let gout = c.geom.maximum_inscribed_circle(c.tol);
            let us = t.elapsed().as_secs_f64() * 1e6;
            if us < geos {
                geos = us;
            }
            if let Ok(g) = gout {
                if let Ok(len) = g.length() {
                    rg = len.max(0.0);
                }
            }
        }
        rows.push(Row {
            name: c.name.clone(),
            verts: c.verts,
            native_us: native,
            geos_us: geos,
            r_native: rn,
            r_geos: rg,
        });
    }

    let mut ratios: Vec<f64> = rows
        .iter()
        .filter(|r| r.native_us > 0.0 && r.geos_us > 0.0)
        .map(|r| r.geos_us / r.native_us)
        .collect();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let total_n: f64 = rows.iter().map(|r| r.native_us).sum();
    let total_g: f64 = rows.iter().map(|r| r.geos_us).sum();

    println!("=== per-case ratio geos/native ===");
    println!(
        "p10 {:.2}x   p25 {:.2}x   median {:.2}x   p75 {:.2}x   p90 {:.2}x",
        pct(&ratios, 0.10),
        pct(&ratios, 0.25),
        pct(&ratios, 0.50),
        pct(&ratios, 0.75),
        pct(&ratios, 0.90)
    );
    println!(
        "\ntotal native {:.1} ms   total geos {:.1} ms   ratio {:.2}x",
        total_n / 1e3,
        total_g / 1e3,
        total_g / total_n
    );
    let native_wins = rows.iter().filter(|r| r.native_us < r.geos_us).count();
    println!(
        "native faster on {native_wins}/{} cases ({:.0}%)",
        rows.len(),
        100.0 * native_wins as f64 / rows.len() as f64
    );

    let buckets = [
        (0usize, 8usize),
        (8, 17),
        (17, 33),
        (33, 65),
        (65, 129),
        (129, 1_000_000),
    ];
    println!("\n=== by exterior vertex count ===");
    println!(
        "{:>10} {:>6} {:>12} {:>12} {:>10}",
        "verts", "n", "native ms", "geos ms", "ratio"
    );
    for (lo, hi) in buckets {
        let sel: Vec<&Row> = rows
            .iter()
            .filter(|r| r.verts >= lo && r.verts < hi)
            .collect();
        if sel.is_empty() {
            continue;
        }
        let n: f64 = sel.iter().map(|r| r.native_us).sum();
        let g: f64 = sel.iter().map(|r| r.geos_us).sum();
        println!(
            "{:>10} {:>6} {:>12.2} {:>12.2} {:>9.2}x",
            format!("{lo}-{hi}"),
            sel.len(),
            n / 1e3,
            g / 1e3,
            g / n
        );
    }

    let mut diffs: Vec<f64> = rows
        .iter()
        .filter(|r| r.r_native > 0.0 && r.r_geos > 0.0)
        .map(|r| (r.r_native - r.r_geos).abs() / r.r_native.max(r.r_geos))
        .collect();
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("\n=== radius agreement ===");
    println!(
        "median {:.2e}   p90 {:.2e}   max {:.2e}",
        pct(&diffs, 0.5),
        pct(&diffs, 0.9),
        diffs[diffs.len() - 1]
    );

    let mut best: Vec<&Row> = rows.iter().collect();
    best.sort_by(|a, b| {
        (b.geos_us / b.native_us)
            .partial_cmp(&(a.geos_us / a.native_us))
            .unwrap()
    });
    println!("\n=== strongest wins for us ===");
    for r in best.iter().take(5) {
        println!(
            "  {:<20} verts={:<4} native {:>8.1} us   geos {:>9.1} us   {:.1}x",
            r.name,
            r.verts,
            r.native_us,
            r.geos_us,
            r.geos_us / r.native_us
        );
    }

    let mut weak: Vec<&Row> = rows.iter().collect();
    weak.sort_by(|a, b| {
        (a.geos_us / a.native_us)
            .partial_cmp(&(b.geos_us / b.native_us))
            .unwrap()
    });
    println!("\n=== weakest cases ===");
    for r in weak.iter().take(5) {
        println!(
            "  {:<20} verts={:<4} native {:>8.1} us   geos {:>9.1} us   {:.2}x",
            r.name,
            r.verts,
            r.native_us,
            r.geos_us,
            r.geos_us / r.native_us
        );
    }
}
