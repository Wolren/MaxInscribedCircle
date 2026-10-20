//! Load the vendored real-world corpus.
//!
//! `benches/data/realworld_290.geojson` is 290 Polish cadastral parcels
//! (BDOT10k, TERYT 2261), the same set IGE benchmarks against. It is vendored so
//! the corpus travels with the crate and a comparison does not depend on a
//! checkout elsewhere on the machine.

#![allow(dead_code)]

use std::path::PathBuf;

use max_inscribed_circle::geo_types::{Coord, LineString, Polygon};

/// Path to the vendored corpus, resolved from the crate root.
pub fn corpus_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("benches")
        .join("data")
        .join("realworld_290.geojson")
}

/// A parcel: its name and polygon.
pub struct Parcel {
    pub name: String,
    pub polygon: Polygon<f64>,
}

/// Read every polygon from the corpus, in file order.
pub fn load_parcels() -> Vec<Parcel> {
    let text = std::fs::read_to_string(corpus_path()).expect("corpus file");
    let root: serde_json::Value = serde_json::from_str(&text).expect("corpus json");

    let mut out = Vec::new();
    let features = root["features"].as_array().expect("features array");
    for (idx, feat) in features.iter().enumerate() {
        let fid = feat["properties"]["fid"]
            .as_i64()
            .map(|v| v.to_string())
            .unwrap_or_else(|| idx.to_string());
        let name = feat["properties"]["TERYT"]
            .as_str()
            .map(|t| format!("fid{fid}/TERYT{t}"))
            .unwrap_or_else(|| format!("fid{fid}"));

        let geom = &feat["geometry"];
        let kind = geom["type"].as_str().unwrap_or("");
        let polys: Vec<Polygon<f64>> = match kind {
            "Polygon" => rings_to_polygon(&geom["coordinates"]).into_iter().collect(),
            "MultiPolygon" => geom["coordinates"]
                .as_array()
                .map(|a| a.iter().filter_map(rings_to_polygon).collect())
                .unwrap_or_default(),
            _ => Vec::new(),
        };
        for (k, p) in polys.into_iter().enumerate() {
            let label = if polys_len(&feat["geometry"]) > 1 {
                format!("{name}#{k}")
            } else {
                name.clone()
            };
            out.push(Parcel {
                name: label,
                polygon: p,
            });
        }
    }
    out
}

fn polys_len(geom: &serde_json::Value) -> usize {
    if geom["type"].as_str() == Some("MultiPolygon") {
        geom["coordinates"].as_array().map(|a| a.len()).unwrap_or(0)
    } else {
        1
    }
}

/// Build a polygon from a GeoJSON coordinate array: `[[ring], [hole], ...]`.
fn rings_to_polygon(coords: &serde_json::Value) -> Option<Polygon<f64>> {
    let rings = coords.as_array()?;
    let mut iter = rings.iter();
    let exterior = ring_from(iter.next()?)?;
    let interiors: Vec<LineString<f64>> = iter.filter_map(ring_from).collect();
    Some(Polygon::new(exterior, interiors))
}

fn ring_from(ring: &serde_json::Value) -> Option<LineString<f64>> {
    let arr = ring.as_array()?;
    let mut pts: Vec<Coord<f64>> = Vec::with_capacity(arr.len());
    for c in arr {
        let pair = c.as_array()?;
        let x = pair.first()?.as_f64()?;
        let y = pair.get(1)?.as_f64()?;
        pts.push(Coord { x, y });
    }
    if pts.len() < 3 {
        return None;
    }
    Some(LineString::from(pts))
}

/// Vertex count of the exterior ring, used to bucket results.
pub fn exterior_len(p: &Polygon<f64>) -> usize {
    p.exterior().0.len()
}
