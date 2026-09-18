//! Solve the maximum inscribed circle of a square and of a multipolygon.
//!
//! Run with `cargo run --example quickstart`.

use max_inscribed_circle::geo_types::{LineString, MultiPolygon, Polygon};
use max_inscribed_circle::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon, MicOptions,
};

fn square(x0: f64, y0: f64, side: f64) -> Polygon<f64> {
    Polygon::new(
        LineString::from(vec![
            (x0, y0),
            (x0 + side, y0),
            (x0 + side, y0 + side),
            (x0, y0 + side),
            (x0, y0),
        ]),
        vec![],
    )
}

fn main() -> Result<(), max_inscribed_circle::MicError> {
    let opts = MicOptions::default();

    let mic = maximum_inscribed_circle(&square(0.0, 0.0, 4.0), &opts)?;
    println!(
        "square: center ({}, {}) radius {} via {:?}/{:?}",
        mic.center.x(),
        mic.center.y(),
        mic.radius,
        mic.used_engine,
        mic.approach
    );
    assert!((mic.center.x() - 2.0).abs() < 1e-9);
    assert!((mic.center.y() - 2.0).abs() < 1e-9);
    assert!((mic.radius - 2.0).abs() < 1e-9);

    let multi = MultiPolygon::new(vec![square(0.0, 0.0, 2.0), square(10.0, 0.0, 6.0)]);
    let best = maximum_inscribed_circle_multipolygon(&multi, &opts)?;
    println!(
        "multipolygon: component {:?} radius {}",
        best.component_index, best.radius
    );
    assert_eq!(best.component_index, Some(1));
    assert!((best.radius - 3.0).abs() < 1e-9);

    Ok(())
}
