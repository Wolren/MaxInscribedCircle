use geo_types::Polygon;
#[cfg(feature = "geos")]
use geos::{Geom, Geometry};

use super::super::index::NearestBoundaryIndex;
use super::super::input::{HostPolygon, SegmentIndex};
use super::super::{MicApproach, MicError, MicOptions, MicResult, MicUsedEngine, RobustMode};

#[cfg(feature = "geos")]
fn polygon_to_wkt(poly: &Polygon<f64>) -> String {
    let mut out = String::from("POLYGON (");
    let ext = poly.exterior();
    out.push('(');
    for (idx, coord) in ext.0.iter().enumerate() {
        if idx > 0 {
            out.push_str(", ");
        }
        out.push_str(&format!("{} {}", coord.x, coord.y));
    }
    out.push(')');
    for hole in poly.interiors() {
        out.push_str(", ");
        out.push('(');
        for (idx, coord) in hole.0.iter().enumerate() {
            if idx > 0 {
                out.push_str(", ");
            }
            out.push_str(&format!("{} {}", coord.x, coord.y));
        }
        out.push(')');
    }
    out.push(')');
    out
}

#[cfg(feature = "geos")]
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn solve_with_geos(
    host: &HostPolygon,
    opts: &MicOptions,
    existing_seg_index: Option<&SegmentIndex>,
) -> std::result::Result<MicResult, MicError> {
    // Use ORIGINAL polygon (not normalized) for GEOS to avoid the bug where
    // normalization breaks certain polygons and causes wrong MIC results
    let poly_wkt = polygon_to_wkt(&host.original_polygon);
    let geom = Geometry::new_from_wkt(&poly_wkt)
        .map_err(|err| MicError::GeosFailed(format!("failed to build GEOS polygon: {err}")))?;

    let tolerance = geos_tolerance(host, opts);
    let radius_geom = geom
        .maximum_inscribed_circle(tolerance)
        .map_err(|err| MicError::GeosFailed(format!("maximum_inscribed_circle failed: {err}")))?;
    let out_wkt = radius_geom
        .to_wkt()
        .map_err(|err| MicError::GeosFailed(format!("failed to decode GEOS MIC output: {err}")))?;

    let (center, boundary_hint) = parse_geos_output(&out_wkt)?;
    let seg_index = match existing_seg_index {
        Some(idx) => idx.clone(),
        None => SegmentIndex::from_host(host),
    };
    let nb_index = NearestBoundaryIndex::new(seg_index, false);
    let Some((nearest_sq, _)) = nb_index.nearest_distance_sq(center.0, center.1) else {
        return Err(MicError::NoCircleFound);
    };

    let line_sq = boundary_hint.map(|bp| {
        let dx = center.0 - bp.0;
        let dy = center.1 - bp.1;
        dx * dx + dy * dy
    });
    let radius_sq = line_sq.unwrap_or(nearest_sq).min(nearest_sq).max(0.0);
    let support_eps = radius_sq.max(1.0) * 1e-10;
    let support_segments = nb_index.supporting_segments(center.0, center.1, radius_sq, support_eps);

    Ok(MicResult {
        center: geo_types::Point::new(center.0, center.1),
        radius: radius_sq.sqrt(),
        radius_sq,
        support_segments,
        candidate_count: 1,
        used_engine: MicUsedEngine::GeosFallback,
        approach: MicApproach::GeosFallback,
        component_index: None,
    })
}

#[cfg(feature = "geos")]
fn geos_tolerance(host: &HostPolygon, opts: &MicOptions) -> f64 {
    let Some((min_x, min_y, max_x, max_y)) = host.bounds() else {
        return 1e-6;
    };
    let diag = (max_x - min_x).hypot(max_y - min_y).max(1.0);
    let factor = match opts.robust_mode {
        RobustMode::FastF64 => 1e-6,
        RobustMode::Filtered => 1e-8,
    };
    (diag * factor).max(1e-12)
}

#[allow(dead_code)]
fn host_polygon_to_wkt(host: &HostPolygon) -> String {
    let mut out = String::from("POLYGON (");
    for ring_id in 0..host.ring_count() {
        if ring_id > 0 {
            out.push_str(", ");
        }
        out.push('(');
        let ring = host.ring_coords(ring_id);
        for (idx, p) in ring.iter().enumerate() {
            if idx > 0 {
                out.push_str(", ");
            }
            out.push_str(&format!("{} {}", p[0], p[1]));
        }
        out.push(')');
    }
    out.push(')');
    out
}

fn parse_geos_output(wkt: &str) -> Result<((f64, f64), Option<(f64, f64)>), MicError> {
    let trimmed = wkt.trim();
    if trimmed.starts_with("LINESTRING") {
        let coords = parse_wkt_coord_list(trimmed)?;
        if coords.len() < 2 {
            return Err(MicError::UnsupportedGeosOutput(format!(
                "LINESTRING has fewer than two coordinates: {trimmed}"
            )));
        }
        return Ok((coords[0], Some(coords[1])));
    }
    if trimmed.starts_with("POINT") {
        let coords = parse_wkt_coord_list(trimmed)?;
        if coords.is_empty() {
            return Err(MicError::UnsupportedGeosOutput(format!(
                "POINT has no coordinates: {trimmed}"
            )));
        }
        return Ok((coords[0], None));
    }

    Err(MicError::UnsupportedGeosOutput(format!(
        "expected LINESTRING/POINT, got: {trimmed}"
    )))
}

fn parse_wkt_coord_list(wkt: &str) -> Result<Vec<(f64, f64)>, MicError> {
    let start = wkt.find('(').ok_or_else(|| {
        MicError::UnsupportedGeosOutput(format!("WKT has no opening parenthesis: {wkt}"))
    })?;
    let end = wkt.rfind(')').ok_or_else(|| {
        MicError::UnsupportedGeosOutput(format!("WKT has no closing parenthesis: {wkt}"))
    })?;
    if end <= start {
        return Err(MicError::UnsupportedGeosOutput(format!(
            "malformed WKT coordinate body: {wkt}"
        )));
    }
    let body = &wkt[start + 1..end];
    if body.trim().is_empty() {
        return Ok(Vec::new());
    }

    let mut out = Vec::new();
    for token in body.split(',') {
        let mut parts = token.split_whitespace();
        let x = parts
            .next()
            .ok_or_else(|| {
                MicError::UnsupportedGeosOutput(format!("missing x coordinate in: {wkt}"))
            })?
            .parse::<f64>()
            .map_err(|err| {
                MicError::UnsupportedGeosOutput(format!("invalid x coordinate: {err}"))
            })?;
        let y = parts
            .next()
            .ok_or_else(|| {
                MicError::UnsupportedGeosOutput(format!("missing y coordinate in: {wkt}"))
            })?
            .parse::<f64>()
            .map_err(|err| {
                MicError::UnsupportedGeosOutput(format!("invalid y coordinate: {err}"))
            })?;
        out.push((x, y));
    }
    Ok(out)
}
