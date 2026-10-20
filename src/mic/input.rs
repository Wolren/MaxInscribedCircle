//! Polygon input normalization and boundary segment tables for the MIC
//! solvers.
//!
//! [`HostPolygon`] flattens `geo_types` rings into one coordinate buffer,
//! and [`SegmentIndex`] exposes the boundary edges as a struct-of-arrays
//! table shared by the indexes and the solvers.
use geo::Area;
use geo_types::{Coord, LineString, Polygon};

use super::MicError;

const NORMALIZE_EPS: f64 = 1e-12;

/// Metadata describing a ring inside the flat coordinate buffer.
#[derive(Debug, Clone)]
pub struct RingMeta {
    /// Index into `HostPolygon::coords` where this ring starts.
    pub start: usize,
    /// Exclusive index into `HostPolygon::coords` where this ring ends.
    pub end: usize,
    /// True for an interior ring, false for the exterior ring.
    pub is_hole: bool,
}

/// Normalized polygon input used by MIC solvers.
#[derive(Debug, Clone)]
pub struct HostPolygon {
    /// Flat coordinate storage for all rings.
    pub coords: Vec<[f64; 2]>,
    /// Ring offsets into `coords`.
    pub rings: Vec<RingMeta>,
    /// The caller's polygon, kept so a comparison solver can hand GEOS the
    /// exact coordinates it was given. Coordinates are never normalized, so
    /// this is the same geometry `coords` holds, and one clone serves both.
    #[cfg(feature = "bench-geos")]
    pub original_polygon: Polygon<f64>,
}

impl HostPolygon {
    /// Flatten a `geo_types` polygon into rings, rejecting non-finite
    /// coordinates and zero-area input.
    pub fn from_polygon(poly: &Polygon<f64>) -> Result<Self, MicError> {
        // Coordinates are used as given, never normalized, so a caller's
        // geometry is solved exactly as supplied. Non-finite values are
        // rejected here because they poison every distance comparison: a NaN
        // never compares greater, so the search keeps an incumbent of zero and
        // returns a circle that has nothing to do with the input.
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

        if coords
            .iter()
            .any(|c| !c[0].is_finite() || !c[1].is_finite())
        {
            return Err(MicError::InvalidInput(
                "polygon has a NaN or infinite coordinate".to_string(),
            ));
        }

        if poly.unsigned_area() <= NORMALIZE_EPS {
            return Err(MicError::InvalidInput("polygon area is zero".to_string()));
        }

        Ok(Self {
            coords,
            rings,
            #[cfg(feature = "bench-geos")]
            original_polygon: poly.clone(),
        })
    }

    /// Coordinates of ring `ring_id`, including its closing point.
    pub fn ring_coords(&self, ring_id: usize) -> &[[f64; 2]] {
        let meta = &self.rings[ring_id];
        &self.coords[meta.start..meta.end]
    }

    /// Coordinates of the exterior ring, including its closing point.
    pub fn outer_ring(&self) -> &[[f64; 2]] {
        self.ring_coords(0)
    }

    /// Number of rings, the exterior ring first.
    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    /// All ring vertices in ring order, with closing points omitted.
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

    /// Bounds as (min_x, min_y, max_x, max_y), or None when there are no coords.
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
    /// Segment start x; index i matches index i of every other array here.
    pub ax: Vec<f64>,
    /// Segment start y, aligned with `ax`.
    pub ay: Vec<f64>,
    /// Segment end x, aligned with `ax`.
    pub bx: Vec<f64>,
    /// Segment end y, aligned with `ax`.
    pub by: Vec<f64>,
    /// Index into `HostPolygon::rings` of the ring owning each segment.
    pub ring_id: Vec<usize>,
    /// Edge position of each segment within its own ring.
    pub edge_id: Vec<usize>,
    /// True when the ring owning the segment is a hole.
    pub is_hole_edge: Vec<bool>,
    /// Minimal x of each segment's bounding box.
    pub bbox_minx: Vec<f64>,
    /// Maximal x of each segment's bounding box.
    pub bbox_maxx: Vec<f64>,
    /// Minimal y of each segment's bounding box.
    pub bbox_miny: Vec<f64>,
    /// Maximal y of each segment's bounding box.
    pub bbox_maxy: Vec<f64>,
    /// Segment direction x, computed as bx - ax.
    pub dir_x: Vec<f64>,
    /// Segment direction y, computed as by - ay.
    pub dir_y: Vec<f64>,
    /// Squared segment length, dir_x * dir_x + dir_y * dir_y.
    pub inv_len_sq: Vec<f64>,
}

impl SegmentIndex {
    /// Build the table from every non-degenerate edge of a host polygon.
    pub fn from_host(host: &HostPolygon) -> Self {
        // Two passes, no intermediates. The first counts the segments that
        // survive the degenerate-edge filter so the columns can be allocated
        // exactly once; the second fills them in place. Building a per-ring Vec
        // of 14-tuples and then flattening it cost three allocations and a
        // repack per solve for a number that is cheap to recount.
        let mut len = 0usize;
        for (rid, meta) in host.rings.iter().enumerate() {
            let ring = host.ring_coords(rid);
            if ring.len() < 2 {
                continue;
            }
            let _ = meta;
            for eid in 0..ring.len() - 1 {
                let a = ring[eid];
                let b = ring[eid + 1];
                let dx = b[0] - a[0];
                let dy = b[1] - a[1];
                if dx * dx + dy * dy > NORMALIZE_EPS {
                    len += 1;
                }
            }
        }

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
            inv_len_sq: Vec::with_capacity(len),
        };

        for (rid, meta) in host.rings.iter().enumerate() {
            let ring = host.ring_coords(rid);
            if ring.len() < 2 {
                continue;
            }
            for eid in 0..ring.len() - 1 {
                let a = ring[eid];
                let b = ring[eid + 1];
                let dx = b[0] - a[0];
                let dy = b[1] - a[1];
                let len_sq = dx * dx + dy * dy;
                if len_sq <= NORMALIZE_EPS {
                    continue;
                }
                index.ax.push(a[0]);
                index.ay.push(a[1]);
                index.bx.push(b[0]);
                index.by.push(b[1]);
                index.ring_id.push(rid);
                index.edge_id.push(eid);
                index.is_hole_edge.push(meta.is_hole);
                index.bbox_minx.push(a[0].min(b[0]));
                index.bbox_maxx.push(a[0].max(b[0]));
                index.bbox_miny.push(a[1].min(b[1]));
                index.bbox_maxy.push(a[1].max(b[1]));
                index.dir_x.push(dx);
                index.dir_y.push(dy);
                index.inv_len_sq.push(1.0 / len_sq);
            }
        }

        index
    }

    /// Number of segments in the table.
    pub fn len(&self) -> usize {
        self.ax.len()
    }

    /// True when the table holds no segments.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Midpoint of segment `seg_idx` as (x, y).
    pub fn midpoint(&self, seg_idx: usize) -> (f64, f64) {
        (
            (self.ax[seg_idx] + self.bx[seg_idx]) * 0.5,
            (self.ay[seg_idx] + self.by[seg_idx]) * 0.5,
        )
    }

    #[inline]
    /// Squared distance from point (x, y) to segment `seg_idx`.
    pub fn point_segment_distance_sq(&self, seg_idx: usize, x: f64, y: f64) -> f64 {
        let ax = self.ax[seg_idx];
        let ay = self.ay[seg_idx];
        let dx = self.dir_x[seg_idx];
        let dy = self.dir_y[seg_idx];
        let inv_len_sq = self.inv_len_sq[seg_idx];

        let t = (((x - ax) * dx + (y - ay) * dy) * inv_len_sq).clamp(0.0, 1.0);
        let px = ax + t * dx;
        let py = ay + t * dy;
        let ex = x - px;
        let ey = y - py;
        ex * ex + ey * ey
    }

    #[cfg(feature = "simd")]
    #[inline]
    /// Smallest squared distance from (x, y) to segments `start..end`.
    pub fn batch_point_segment_distance_sq(
        &self,
        px: f64,
        py: f64,
        start: usize,
        end: usize,
    ) -> f64 {
        use wide::f64x4;

        let x = f64x4::splat(px);
        let y = f64x4::splat(py);

        let mut best = f64::INFINITY;

        let mut i = start;
        while i + 4 <= end {
            let ax = f64x4::from([self.ax[i], self.ax[i + 1], self.ax[i + 2], self.ax[i + 3]]);
            let ay = f64x4::from([self.ay[i], self.ay[i + 1], self.ay[i + 2], self.ay[i + 3]]);
            let dx = f64x4::from([
                self.dir_x[i],
                self.dir_x[i + 1],
                self.dir_x[i + 2],
                self.dir_x[i + 3],
            ]);
            let dy = f64x4::from([
                self.dir_y[i],
                self.dir_y[i + 1],
                self.dir_y[i + 2],
                self.dir_y[i + 3],
            ]);
            let ls = f64x4::from([
                self.inv_len_sq[i],
                self.inv_len_sq[i + 1],
                self.inv_len_sq[i + 2],
                self.inv_len_sq[i + 3],
            ]);

            let t = ((x - ax) * dx + (y - ay) * dy) * ls;
            let t = t.max(f64x4::ZERO).min(f64x4::ONE);
            let pdx = ax + t * dx;
            let pdy = ay + t * dy;
            let edx = x - pdx;
            let edy = y - pdy;
            let dist_sq = edx * edx + edy * edy;

            best = best.min(dist_sq.to_array().into_iter().fold(f64::INFINITY, f64::min));
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
    /// Smallest squared distance from (x, y) to segments `start..end`.
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
    /// Smallest squared distance over `start..end` and the segment that achieves it.
    pub fn batch_point_segment_distance_sq_with_index(
        &self,
        px: f64,
        py: f64,
        start: usize,
        end: usize,
    ) -> (f64, usize) {
        use wide::f64x4;

        let x = f64x4::splat(px);
        let y = f64x4::splat(py);

        let mut best_dist = f64::INFINITY;
        let mut best_idx = start;

        let mut i = start;
        while i + 4 <= end {
            let ax = f64x4::from([self.ax[i], self.ax[i + 1], self.ax[i + 2], self.ax[i + 3]]);
            let ay = f64x4::from([self.ay[i], self.ay[i + 1], self.ay[i + 2], self.ay[i + 3]]);
            let dx = f64x4::from([
                self.dir_x[i],
                self.dir_x[i + 1],
                self.dir_x[i + 2],
                self.dir_x[i + 3],
            ]);
            let dy = f64x4::from([
                self.dir_y[i],
                self.dir_y[i + 1],
                self.dir_y[i + 2],
                self.dir_y[i + 3],
            ]);
            let ls = f64x4::from([
                self.inv_len_sq[i],
                self.inv_len_sq[i + 1],
                self.inv_len_sq[i + 2],
                self.inv_len_sq[i + 3],
            ]);

            let t = ((x - ax) * dx + (y - ay) * dy) * ls;
            let t = t.max(f64x4::ZERO).min(f64x4::ONE);
            let pdx = ax + t * dx;
            let pdy = ay + t * dy;
            let edx = x - pdx;
            let edy = y - pdy;
            let dist_sq = edx * edx + edy * edy;
            let lanes = dist_sq.to_array();

            // Scan the four lanes against the incumbent, not against the lane
            // minimum, so the reported index is the one that holds the distance.
            for (lane, d) in lanes.iter().enumerate() {
                if *d < best_dist {
                    best_dist = *d;
                    best_idx = i + lane;
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
    /// Smallest squared distance over `start..end` and the segment that achieves it.
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
