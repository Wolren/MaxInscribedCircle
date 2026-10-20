//! Accuracy of the real-corpus answer against independently measured truth.
//!
//! For every parcel: take the solver's centre, re-measure clearance against the
//! segment table (must equal the reported radius), and compare the radius to
//! the best clearance found by a dense interior scan. The scan is a lower bound
//! on the true MIC, so `radius < scan` would be a definite error and
//! `radius > measured clearance` would be an unsound circle.
//!
//! Run: cargo run --release --example real_accuracy

use max_inscribed_circle::mic::index::{NearestBoundaryIndex, PipIndex};
use max_inscribed_circle::mic::input::{HostPolygon, SegmentIndex};
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

#[path = "../benches/corpus.rs"]
mod corpus;

fn main() {
    let parcels = corpus::load_parcels();
    let mut worst_short = 0.0f64;
    let mut worst_short_name = String::new();
    let mut worst_over = 0.0f64;
    let mut over_names = Vec::new();
    let mut total_short = 0.0f64;
    let mut n = 0usize;

    for p in &parcels {
        let Ok(r) = maximum_inscribed_circle(&p.polygon, &MicOptions::default()) else {
            continue;
        };
        let Ok(host) = HostPolygon::from_polygon(&p.polygon) else {
            continue;
        };
        let seg = SegmentIndex::from_host(&host);
        let pip = PipIndex::new(&host);
        let nb = NearestBoundaryIndex::new(seg);

        // 1. the reported radius must equal the clearance at its own centre.
        let measured = nb
            .nearest_distance_sq(r.center.x(), r.center.y())
            .map(|(d, _)| d.max(0.0).sqrt())
            .unwrap_or(0.0);
        let over = (r.radius - measured) / measured.max(1e-300);
        if over > 1e-12 {
            if over > worst_over {
                worst_over = over;
            }
            over_names.push(p.name.clone());
        }

        // 2. how far short of the best sampled clearance the answer falls.
        let (min_x, min_y, max_x, max_y) = host.bounds().unwrap();
        let steps = 220;
        let mut best = 0.0f64;
        for i in 0..=steps {
            let x = min_x + (max_x - min_x) * (i as f64) / steps as f64;
            for j in 0..=steps {
                let y = min_y + (max_y - min_y) * (j as f64) / steps as f64;
                if !pip.contains_strict_xy(x, y) {
                    continue;
                }
                if let Some((d, _)) = nb.nearest_distance_sq(x, y) {
                    let d = d.max(0.0).sqrt();
                    if d > best {
                        best = d;
                    }
                }
            }
        }
        let short = (best - r.radius) / best.max(1e-300);
        total_short += short.max(0.0);
        if short > worst_short {
            worst_short = short;
            worst_short_name = p.name.clone();
        }
        n += 1;
    }

    println!("parcels checked: {n}");
    println!(
        "radius above its own measured clearance: {} (worst {:+.3e})",
        over_names.len(),
        worst_over
    );
    if !over_names.is_empty() {
        println!("  offenders: {:?}", &over_names[..over_names.len().min(6)]);
    }
    println!(
        "worst shortfall vs dense scan: {:.3e} relative ({})",
        worst_short, worst_short_name
    );
    println!(
        "mean shortfall vs dense scan: {:.3e} relative",
        total_short / n as f64
    );
}
