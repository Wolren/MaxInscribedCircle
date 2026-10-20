//! georust ecosystem interop (feature `geo-traits`).
//!
//! Two surfaces. First, the solvers accept any foreign polygon that
//! implements `geo_traits::PolygonTrait` (wkb views, geoarrow arrays,
//! geozero streams) without materializing a `geo_types::Polygon` first;
//! conversion uses `to_polygon()` and keeps only the first two dimensions.
//! Second, [`maximum_inscribed_circle_geometry`] dispatches over
//! `geo_types::Geometry` so callers holding mixed geometry enums get one
//! entry point.
//!
//! MultiPolygon policy: the native multipolygon solve runs and keeps the
//! best component, exactly as the Rust API does. Other geometry variants
//! are rejected with an error naming the variant, never silently
//! mis-solved.

use super::mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon, MicError, MicOptions,
    MicResult,
};
use geo_traits::to_geo::ToGeoPolygon;
use geo_traits::PolygonTrait;
use geo_types::{Geometry, MultiPolygon, Polygon};

/// Solve MIC on any foreign polygon implementing `PolygonTrait`.
pub fn maximum_inscribed_circle_polygon_trait<G>(
    geom: &G,
    opts: &MicOptions,
) -> Result<MicResult, MicError>
where
    G: PolygonTrait<T = f64>,
{
    let poly: Polygon<f64> = geom.to_polygon();
    maximum_inscribed_circle(&poly, opts)
}

/// Solve MIC on any foreign multipolygon implementing `MultiPolygonTrait`.
pub fn maximum_inscribed_circle_multipolygon_trait<G>(
    geom: &G,
    opts: &MicOptions,
) -> Result<MicResult, MicError>
where
    G: geo_traits::MultiPolygonTrait<T = f64>,
{
    use geo_traits::to_geo::ToGeoMultiPolygon;
    let multi: MultiPolygon<f64> = geom.to_multi_polygon();
    maximum_inscribed_circle_multipolygon(&multi, opts)
}

fn variant_name(geom: &Geometry<f64>) -> &'static str {
    match geom {
        Geometry::Point(_) => "Point",
        Geometry::Line(_) => "Line",
        Geometry::LineString(_) => "LineString",
        Geometry::Polygon(_) => "Polygon",
        Geometry::MultiPoint(_) => "MultiPoint",
        Geometry::MultiLineString(_) => "MultiLineString",
        Geometry::MultiPolygon(_) => "MultiPolygon",
        Geometry::GeometryCollection(_) => "GeometryCollection",
        Geometry::Rect(_) => "Rect",
        Geometry::Triangle(_) => "Triangle",
    }
}

/// Solve MIC on a `geo_types::Geometry`.
///
/// Polygon solves directly. MultiPolygon runs the native multipolygon
/// solve and keeps the best component. Every other variant is rejected
/// with an error naming it.
pub fn maximum_inscribed_circle_geometry(
    geom: &Geometry<f64>,
    opts: &MicOptions,
) -> Result<MicResult, MicError> {
    match geom {
        Geometry::Polygon(p) => maximum_inscribed_circle(p, opts),
        Geometry::MultiPolygon(mp) => maximum_inscribed_circle_multipolygon(mp, opts),
        other => Err(MicError::InvalidInput(format!(
            "maximum_inscribed_circle_geometry needs Polygon or MultiPolygon, got {}",
            variant_name(other)
        ))),
    }
}
