//! geo-traits interop tests: foreign polygons solve identically to native.
#![cfg(feature = "geo-traits")]

use geo_types::{Coord, LineString, MultiPolygon, Polygon};
use max_inscribed_circle::geo_interop::{
    maximum_inscribed_circle_geometry, maximum_inscribed_circle_polygon_trait,
};
use max_inscribed_circle::{maximum_inscribed_circle, MicOptions};

fn square() -> Polygon<f64> {
    Polygon::new(
        LineString::from(vec![
            Coord { x: 0.0, y: 0.0 },
            Coord { x: 10.0, y: 0.0 },
            Coord { x: 10.0, y: 10.0 },
            Coord { x: 0.0, y: 10.0 },
            Coord { x: 0.0, y: 0.0 },
        ]),
        vec![],
    )
}

#[test]
fn polygon_trait_matches_native() {
    let poly = square();
    let opts = MicOptions::default();
    let native = maximum_inscribed_circle(&poly, &opts).unwrap();
    // `&Polygon` implements `PolygonTrait` via geo-traits' geo-types impls,
    // exercising the foreign-type path with a borrowed view.
    let foreign = maximum_inscribed_circle_polygon_trait(&&poly, &opts).unwrap();
    assert!((native.radius - foreign.radius).abs() < 1e-9);
}

#[test]
fn geometry_dispatch_polygon_and_multi() {
    let opts = MicOptions::default();
    let geom = geo_types::Geometry::Polygon(square());
    let r = maximum_inscribed_circle_geometry(&geom, &opts).unwrap();
    assert!((r.radius - 5.0).abs() < 1e-6);

    let small = Polygon::new(
        LineString::from(vec![
            Coord { x: 20.0, y: 20.0 },
            Coord { x: 22.0, y: 20.0 },
            Coord { x: 22.0, y: 22.0 },
            Coord { x: 20.0, y: 22.0 },
            Coord { x: 20.0, y: 20.0 },
        ]),
        vec![],
    );
    let geom = geo_types::Geometry::MultiPolygon(MultiPolygon(vec![small, square()]));
    let r = maximum_inscribed_circle_geometry(&geom, &opts).unwrap();
    assert!((r.radius - 5.0).abs() < 1e-6);
    assert_eq!(r.component_index, Some(1));
}

#[test]
fn geometry_dispatch_rejects_point() {
    let geom = geo_types::Geometry::Point(geo_types::Point::new(1.0, 2.0));
    let err =
        maximum_inscribed_circle_geometry(&geom, &MicOptions::default()).unwrap_err();
    assert!(err.to_string().contains("Point"));
}
