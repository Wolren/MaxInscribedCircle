//! Point-in-polygon index over the rings of a host polygon.
//!
//! Containment uses even-odd ray casting. Ring bounding boxes short-circuit the
//! rings that cannot contain the query, which matters on polygons with holes:
//! the point only pays for the rings whose box it falls inside.
use super::super::input::HostPolygon;

/// Point-in-polygon test with hole support.
#[derive(Debug, Clone)]
pub struct PipIndex {
    coords: Vec<[f64; 2]>,
    rings: Vec<RingMeta>,
    ring_bboxes: Vec<RingBbox>,
}

#[derive(Debug, Clone)]
struct RingMeta {
    start: usize,
    end: usize,
    is_hole: bool,
}

#[derive(Debug, Clone)]
struct RingBbox {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
}

impl PipIndex {
    /// Build the containment index from the rings of a host polygon.
    pub fn new(host: &HostPolygon) -> Self {
        let mut rings = Vec::with_capacity(host.rings.len());
        let mut ring_bboxes = Vec::with_capacity(host.rings.len());

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

    #[inline]
    /// True when (x, y) is inside the polygon, with hole interiors excluded.
    pub fn contains_strict_xy(&self, x: f64, y: f64) -> bool {
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

fn point_in_ring(x: f64, y: f64, ring: &[[f64; 2]]) -> bool {
    #[cfg(not(feature = "shewchuk"))]
    {
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
    #[cfg(feature = "shewchuk")]
    {
        use robust::{orient2d, Coord};
        let q = Coord { x, y };
        let mut inside = false;
        for w in ring.windows(2) {
            let a = Coord {
                x: w[0][0],
                y: w[0][1],
            };
            let b = Coord {
                x: w[1][0],
                y: w[1][1],
            };
            // Same half-open rule as the ray cast: a vertex on the ray counts once.
            if (a.y > y) != (b.y > y) {
                let o = orient2d(a, b, q);
                let left = if b.y > a.y { o > 0.0 } else { o < 0.0 };
                if left {
                    inside = !inside;
                }
            }
        }
        inside
    }
}
