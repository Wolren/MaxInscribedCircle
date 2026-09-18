//! Point-in-polygon index: exact winding numbers behind the `shewchuk`
//! feature, even-odd ray casting otherwise.
#[cfg(feature = "shewchuk")]
use crate::winding_index::WindingIndex;
#[cfg(feature = "shewchuk")]
use geo_types::{Coord, LineString, Polygon};

/// Point-in-polygon test.
///
/// Without the `shewchuk` feature: classic even-odd ray casting (fast, but may
/// have edge cases near vertices).
///
/// With the `shewchuk` feature: exact winding number using Shewchuk adaptive
/// arithmetic + x-interval spatial acceleration.
#[derive(Debug, Clone)]
pub struct PipIndex {
    #[cfg(feature = "shewchuk")]
    winding: WindingIndex,
    /// Flat vertex data for ray-casting fallback.
    #[cfg(not(feature = "shewchuk"))]
    coords: Vec<[f64; 2]>,
    #[cfg(not(feature = "shewchuk"))]
    rings: Vec<RingMeta>,
    #[cfg(not(feature = "shewchuk"))]
    ring_bboxes: Vec<RingBbox>,
}

#[cfg(not(feature = "shewchuk"))]
#[derive(Debug, Clone)]
struct RingMeta {
    start: usize,
    end: usize,
    is_hole: bool,
}

#[cfg(not(feature = "shewchuk"))]
#[derive(Debug, Clone)]
struct RingBbox {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
}

impl PipIndex {
    /// Build the containment index from the rings of a host polygon.
    pub fn new(host: &super::super::input::HostPolygon) -> Self {
        #[cfg(feature = "shewchuk")]
        {
            let exterior = LineString::from(
                host.coords
                    .iter()
                    .take(host.rings[0].end)
                    .map(|c| Coord { x: c[0], y: c[1] })
                    .collect::<Vec<_>>(),
            );

            let mut interiors = Vec::new();
            for ring_meta in host.rings.iter().skip(1) {
                let ring_coords: Vec<Coord<f64>> = host.coords[ring_meta.start..ring_meta.end]
                    .iter()
                    .map(|c| Coord { x: c[0], y: c[1] })
                    .collect();
                interiors.push(LineString::from(ring_coords));
            }

            let poly = Polygon::new(exterior, interiors);
            return Self {
                winding: WindingIndex::from_polygon(&poly),
            };
        }
        #[cfg(not(feature = "shewchuk"))]
        {
            let mut rings = Vec::new();
            let mut ring_bboxes = Vec::new();

            for ring_meta in &host.rings {
                let coords = &host.coords[ring_meta.start..ring_meta.end];
                let mut x_min = f64::INFINITY;
                let mut x_max = f64::NEG_INFINITY;
                let mut y_min = f64::INFINITY;
                let mut y_max = f64::NEG_INFINITY;
                for c in coords {
                    x_min = x_min.min(c[0]);
                    x_max = x_max.max(c[0]);
                    y_min = y_min.min(c[1]);
                    y_max = y_max.max(c[1]);
                }
                ring_bboxes.push(RingBbox {
                    x_min,
                    x_max,
                    y_min,
                    y_max,
                });
                rings.push(RingMeta {
                    start: ring_meta.start,
                    end: ring_meta.end,
                    is_hole: ring_meta.is_hole,
                });
            }

            Self {
                coords: host.coords.clone(),
                rings,
                ring_bboxes,
            }
        }
    }

    #[inline]
    /// True when (x, y) is inside the polygon, with hole interiors excluded.
    pub fn contains_strict_xy(&self, x: f64, y: f64) -> bool {
        #[cfg(feature = "shewchuk")]
        {
            return self.winding.contains(x, y);
        }
        #[cfg(not(feature = "shewchuk"))]
        {
            // Even-odd ray casting with ring-bbox pre-filter
            let mut inside = false;
            for (ring_meta, bbox) in self.rings.iter().zip(self.ring_bboxes.iter()) {
                if x < bbox.x_min || x > bbox.x_max || y < bbox.y_min || y > bbox.y_max {
                    continue;
                }
                let ring_in = point_in_ring(x, y, &self.coords[ring_meta.start..ring_meta.end]);
                if ring_meta.is_hole {
                    if ring_in {
                        inside = false;
                    }
                } else {
                    inside ^= ring_in;
                }
            }
            inside
        }
    }
}

#[cfg(not(feature = "shewchuk"))]
fn point_in_ring(x: f64, y: f64, ring: &[[f64; 2]]) -> bool {
    let mut inside = false;
    for w in ring.windows(2) {
        let (ax, ay) = (w[0][0], w[0][1]);
        let (bx, by) = (w[1][0], w[1][1]);
        if (ay > y) != (by > y) {
            let intersect_x = ax + (y - ay) * (bx - ax) / (by - ay);
            if x < intersect_x {
                inside = !inside;
            }
        }
    }
    inside
}
