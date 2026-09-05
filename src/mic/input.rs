use geo::Area;
use geo_types::{Coord, LineString, Polygon};

use super::MicError;

const NORMALIZE_EPS: f64 = 1e-12;

/// Metadata describing a ring inside the flat coordinate buffer.
#[derive(Debug, Clone)]
pub struct RingMeta {
    pub start: usize,
    pub end: usize,
    pub is_hole: bool,
}

/// Normalized polygon input used by MIC solvers.
#[derive(Debug, Clone)]
pub struct HostPolygon {
    /// Flat coordinate storage for all rings.
    pub coords: Vec<[f64; 2]>,
    /// Ring offsets into `coords`.
    pub rings: Vec<RingMeta>,
    /// Canonicalized geometry used for predicates.
    pub polygon: Polygon<f64>,
    /// Original polygon before normalization (for GEOS to get correct results)
    pub original_polygon: Polygon<f64>,
}

impl HostPolygon {
    pub fn from_polygon(poly: &Polygon<f64>) -> Result<Self, MicError> {
        // NO normalization - use original polygon coordinates directly to avoid any geometry alteration
        let mut coords = Vec::new();
        let mut rings = Vec::with_capacity(1 + poly.interiors().len());

        // Exterior ring - store with closing point
        let ext_coords: Vec<[f64; 2]> = poly.exterior().0.iter().map(|c| [c.x, c.y]).collect();
        let ext_start = coords.len();
        coords.extend_from_slice(&ext_coords);
        let ext_end = coords.len();
        rings.push(RingMeta {
            start: ext_start,
            end: ext_end,
            is_hole: false,
        });

        // Interior rings (holes)
        for hole in poly.interiors() {
            let hole_coords: Vec<[f64; 2]> = hole.0.iter().map(|c| [c.x, c.y]).collect();
            let start = coords.len();
            coords.extend_from_slice(&hole_coords);
            let end = coords.len();
            rings.push(RingMeta {
                start,
                end,
                is_hole: true,
            });
        }

        if poly.unsigned_area() <= NORMALIZE_EPS {
            return Err(MicError::InvalidInput("polygon area is zero".to_string()));
        }

        Ok(Self {
            coords,
            rings,
            polygon: poly.clone(),
            original_polygon: poly.clone(),
        })
    }

    pub fn ring_coords(&self, ring_id: usize) -> &[[f64; 2]] {
        let meta = &self.rings[ring_id];
        &self.coords[meta.start..meta.end]
    }

    pub fn outer_ring(&self) -> &[[f64; 2]] {
        self.ring_coords(0)
    }

    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    pub fn unique_vertices(&self) -> Vec<[f64; 2]> {
        let mut out = Vec::new();
        for ring_id in 0..self.rings.len() {
            let ring = self.ring_coords(ring_id);
            if ring.len() < 2 {
                continue;
            }
            for p in &ring[..ring.len() - 1] {
                out.push(*p);
            }
        }
        out
    }

    pub fn bounds(&self) -> Option<(f64, f64, f64, f64)> {
        if self.coords.is_empty() {
            return None;
        }
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;

        for p in &self.coords {
            min_x = min_x.min(p[0]);
            min_y = min_y.min(p[1]);
            max_x = max_x.max(p[0]);
            max_y = max_y.max(p[1]);
        }

        Some((min_x, min_y, max_x, max_y))
    }
}

/// Struct-of-arrays segment table over all ring edges.
#[derive(Debug, Clone)]
pub struct SegmentIndex {
    pub ax: Vec<f64>,
    pub ay: Vec<f64>,
    pub bx: Vec<f64>,
    pub by: Vec<f64>,
    pub ring_id: Vec<usize>,
    pub edge_id: Vec<usize>,
    pub is_hole_edge: Vec<bool>,
    pub bbox_minx: Vec<f64>,
    pub bbox_maxx: Vec<f64>,
    pub bbox_miny: Vec<f64>,
    pub bbox_maxy: Vec<f64>,
    pub dir_x: Vec<f64>,
    pub dir_y: Vec<f64>,
    pub len_sq: Vec<f64>,
}

impl SegmentIndex {
    pub fn from_host(host: &HostPolygon) -> Self {
        // Collect ring data and flatten sequentially
        let ring_data: Vec<_> = host
            .rings
            .iter()
            .enumerate()
            .map(|(rid, meta)| {
                let ring = host.ring_coords(rid);
                let mut edges = Vec::new();
                if ring.len() >= 2 {
                    for eid in 0..ring.len() - 1 {
                        let a = ring[eid];
                        let b = ring[eid + 1];
                        let dx = b[0] - a[0];
                        let dy = b[1] - a[1];
                        let len_sq = dx * dx + dy * dy;
                        if len_sq <= NORMALIZE_EPS {
                            continue;
                        }
                        edges.push((
                            a[0],
                            a[1],
                            b[0],
                            b[1],
                            rid,
                            eid,
                            meta.is_hole,
                            a[0].min(b[0]),
                            a[0].max(b[0]),
                            a[1].min(b[1]),
                            a[1].max(b[1]),
                            dx,
                            dy,
                            len_sq,
                        ));
                    }
                }
                edges
            })
            .collect();

        let flat_results: Vec<_> = ring_data.iter().flat_map(|edges| edges.iter()).collect();

        // Allocate exact capacity and pack
        let len = flat_results.len();
        let mut index = Self {
            ax: Vec::with_capacity(len),
            ay: Vec::with_capacity(len),
            bx: Vec::with_capacity(len),
            by: Vec::with_capacity(len),
            ring_id: Vec::with_capacity(len),
            edge_id: Vec::with_capacity(len),
            is_hole_edge: Vec::with_capacity(len),
            bbox_minx: Vec::with_capacity(len),
            bbox_maxx: Vec::with_capacity(len),
            bbox_miny: Vec::with_capacity(len),
            bbox_maxy: Vec::with_capacity(len),
            dir_x: Vec::with_capacity(len),
            dir_y: Vec::with_capacity(len),
            len_sq: Vec::with_capacity(len),
        };

        for tup in flat_results {
            index.ax.push(tup.0);
            index.ay.push(tup.1);
            index.bx.push(tup.2);
            index.by.push(tup.3);
            index.ring_id.push(tup.4);
            index.edge_id.push(tup.5);
            index.is_hole_edge.push(tup.6);
            index.bbox_minx.push(tup.7);
            index.bbox_maxx.push(tup.8);
            index.bbox_miny.push(tup.9);
            index.bbox_maxy.push(tup.10);
            index.dir_x.push(tup.11);
            index.dir_y.push(tup.12);
            index.len_sq.push(tup.13);
        }

        index
    }

    pub fn len(&self) -> usize {
        self.ax.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn midpoint(&self, seg_idx: usize) -> (f64, f64) {
        (
            (self.ax[seg_idx] + self.bx[seg_idx]) * 0.5,
            (self.ay[seg_idx] + self.by[seg_idx]) * 0.5,
        )
    }

    #[inline]
    pub fn point_segment_distance_sq(&self, seg_idx: usize, x: f64, y: f64) -> f64 {
        let ax = self.ax[seg_idx];
        let ay = self.ay[seg_idx];
        let dx = self.dir_x[seg_idx];
        let dy = self.dir_y[seg_idx];
        let len_sq = self.len_sq[seg_idx];

        let t = (((x - ax) * dx + (y - ay) * dy) / len_sq).clamp(0.0, 1.0);
        let px = ax + t * dx;
        let py = ay + t * dy;
        let ex = x - px;
        let ey = y - py;
        ex * ex + ey * ey
    }

    #[cfg(feature = "simd")]
    #[inline]
    pub fn batch_point_segment_distance_sq(
        &self,
        px: f64,
        py: f64,
        start: usize,
        end: usize,
    ) -> f64 {
        use std::simd::f64x4;
        use std::simd::num::SimdFloat;

        let x = f64x4::splat(px);
        let y = f64x4::splat(py);

        let mut best = f64::INFINITY;

        let mut i = start;
        while i + 4 <= end {
            let ax =
                f64x4::from_array([self.ax[i], self.ax[i + 1], self.ax[i + 2], self.ax[i + 3]]);
            let ay =
                f64x4::from_array([self.ay[i], self.ay[i + 1], self.ay[i + 2], self.ay[i + 3]]);
            let dx = f64x4::from_array([
                self.dir_x[i],
                self.dir_x[i + 1],
                self.dir_x[i + 2],
                self.dir_x[i + 3],
            ]);
            let dy = f64x4::from_array([
                self.dir_y[i],
                self.dir_y[i + 1],
                self.dir_y[i + 2],
                self.dir_y[i + 3],
            ]);
            let ls = f64x4::from_array([
                self.len_sq[i],
                self.len_sq[i + 1],
                self.len_sq[i + 2],
                self.len_sq[i + 3],
            ]);

            let t = ((x - ax) * dx + (y - ay) * dy) / ls;
            let t = t.simd_clamp(f64x4::splat(0.0), f64x4::splat(1.0));
            let pdx = ax + t * dx;
            let pdy = ay + t * dy;
            let edx = x - pdx;
            let edy = y - pdy;
            let dist_sq = edx * edx + edy * edy;

            best = best.min(dist_sq.reduce_min());
            i += 4;
        }

        for j in i..end {
            let d = self.point_segment_distance_sq(j, px, py);
            if d < best {
                best = d;
            }
        }

        best
    }

    #[cfg(not(feature = "simd"))]
    pub fn batch_point_segment_distance_sq(&self, x: f64, y: f64, start: usize, end: usize) -> f64 {
        let mut best = f64::INFINITY;
        for i in start..end {
            let d = self.point_segment_distance_sq(i, x, y);
            if d < best {
                best = d;
            }
        }
        best
    }

    #[cfg(feature = "simd")]
    #[inline]
    pub fn batch_point_segment_distance_sq_with_index(
        &self,
        px: f64,
        py: f64,
        start: usize,
        end: usize,
    ) -> (f64, usize) {
        use std::simd::f64x4;
        use std::simd::num::SimdFloat;

        let x = f64x4::splat(px);
        let y = f64x4::splat(py);

        let mut best_dist = f64::INFINITY;
        let mut best_idx = start;

        let mut i = start;
        while i + 4 <= end {
            let ax =
                f64x4::from_array([self.ax[i], self.ax[i + 1], self.ax[i + 2], self.ax[i + 3]]);
            let ay =
                f64x4::from_array([self.ay[i], self.ay[i + 1], self.ay[i + 2], self.ay[i + 3]]);
            let dx = f64x4::from_array([
                self.dir_x[i],
                self.dir_x[i + 1],
                self.dir_x[i + 2],
                self.dir_x[i + 3],
            ]);
            let dy = f64x4::from_array([
                self.dir_y[i],
                self.dir_y[i + 1],
                self.dir_y[i + 2],
                self.dir_y[i + 3],
            ]);
            let ls = f64x4::from_array([
                self.len_sq[i],
                self.len_sq[i + 1],
                self.len_sq[i + 2],
                self.len_sq[i + 3],
            ]);

            let t = ((x - ax) * dx + (y - ay) * dy) / ls;
            let t = t.simd_clamp(f64x4::splat(0.0), f64x4::splat(1.0));
            let pdx = ax + t * dx;
            let pdy = ay + t * dy;
            let edx = x - pdx;
            let edy = y - pdy;
            let dist_sq = edx * edx + edy * edy;

            let batch_min = dist_sq.reduce_min();
            if batch_min < best_dist {
                best_dist = batch_min;
                // Find which lane has the minimum
                for lane in 0..4 {
                    let d = dist_sq[lane];
                    if d < best_dist {
                        best_dist = d;
                        best_idx = i + lane;
                    }
                }
            }
            i += 4;
        }

        for j in i..end {
            let d = self.point_segment_distance_sq(j, px, py);
            if d < best_dist {
                best_dist = d;
                best_idx = j;
            }
        }

        (best_dist, best_idx)
    }

    #[cfg(not(feature = "simd"))]
    #[inline]
    pub fn batch_point_segment_distance_sq_with_index(
        &self,
        x: f64,
        y: f64,
        start: usize,
        end: usize,
    ) -> (f64, usize) {
        let mut best = f64::INFINITY;
        let mut best_idx = start;
        for i in start..end {
            let d = self.point_segment_distance_sq(i, x, y);
            if d < best {
                best = d;
                best_idx = i;
            }
        }
        (best, best_idx)
    }
}

#[allow(dead_code)]
fn normalize_ring(ring: &LineString<f64>, is_hole: bool) -> Result<Vec<[f64; 2]>, MicError> {
    let mut pts = Vec::<[f64; 2]>::new();
    for c in &ring.0 {
        if !c.x.is_finite() || !c.y.is_finite() {
            return Err(MicError::InvalidInput(
                "ring contains non-finite coordinates".to_string(),
            ));
        }
        let p = [c.x, c.y];
        if pts
            .last()
            .map(|last| approx_same(*last, p))
            .unwrap_or(false)
        {
            continue;
        }
        pts.push(p);
    }

    if pts.len() < 3 {
        return Err(MicError::InvalidInput(
            "ring has fewer than 3 distinct vertices".to_string(),
        ));
    }

    if approx_same(
        *pts.first().expect("ring has first"),
        *pts.last().expect("ring has last"),
    ) {
        pts.pop();
    }

    if pts.len() < 3 {
        return Err(MicError::InvalidInput(
            "ring collapsed after closure normalization".to_string(),
        ));
    }

    let signed_area = ring_signed_area_open(&pts);
    if signed_area.abs() <= NORMALIZE_EPS {
        return Err(MicError::InvalidInput(
            "ring area is zero after normalization".to_string(),
        ));
    }

    let should_be_ccw = !is_hole;
    let is_ccw = signed_area > 0.0;
    if should_be_ccw != is_ccw {
        pts.reverse();
    }

    pts.push(*pts.first().expect("normalized ring has first"));
    Ok(pts)
}

#[allow(dead_code)]
fn ring_to_linestring(ring: &[[f64; 2]]) -> LineString<f64> {
    let coords = ring
        .iter()
        .map(|p| Coord { x: p[0], y: p[1] })
        .collect::<Vec<_>>();
    LineString::from(coords)
}

#[allow(dead_code)]
fn approx_same(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] - b[0]).abs() <= NORMALIZE_EPS && (a[1] - b[1]).abs() <= NORMALIZE_EPS
}

#[allow(dead_code)]
fn ring_signed_area_open(open_ring: &[[f64; 2]]) -> f64 {
    let n = open_ring.len();
    let mut sum = 0.0;
    for i in 0..n {
        let a = open_ring[i];
        let b = open_ring[(i + 1) % n];
        sum += a[0] * b[1] - b[0] * a[1];
    }
    sum * 0.5
}
