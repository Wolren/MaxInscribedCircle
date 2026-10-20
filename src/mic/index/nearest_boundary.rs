//! Nearest-boundary distance queries over a segment table.
//!
//! Backed by [`Bvh`]. A segment table small enough to scan outright is scanned,
//! because a tree walk costs more than it saves there.
use super::super::input::SegmentIndex;
use super::bvh::Bvh;

/// Below this many segments a linear scan beats a tree descent.
const LINEAR_SCAN_MAX: usize = 32;

/// Nearest-boundary distance queries over a segment table.
#[derive(Debug, Clone)]
pub struct NearestBoundaryIndex {
    segments: SegmentIndex,
    bvh: Bvh,
}

impl NearestBoundaryIndex {
    /// Build the index over `segments`.
    pub fn new(segments: SegmentIndex) -> Self {
        let bvh = Bvh::new(&segments);
        Self { segments, bvh }
    }

    #[inline]
    /// The segment table the index was built from.
    pub fn segments(&self) -> &SegmentIndex {
        &self.segments
    }

    #[inline]
    /// Squared distance to the closest boundary and its segment index, if any.
    pub fn nearest_distance_sq(&self, x: f64, y: f64) -> Option<(f64, usize)> {
        if self.segments.is_empty() {
            return None;
        }
        if self.segments.len() <= LINEAR_SCAN_MAX {
            return Some(linear_scan_nearest(&self.segments, x, y));
        }
        self.bvh.nearest_distance_sq(&self.segments, x, y)
    }

    /// Nearest distance, or `None` once a segment is proven within `threshold_sq`.
    ///
    /// A caller that only needs to know whether a point can beat an incumbent
    /// stops caring about the exact distance as soon as one segment is closer
    /// than the bar. `seed` is a nearby segment index to try first, typically
    /// the support of the parent cell.
    #[inline]
    pub fn nearest_distance_sq_above(
        &self,
        x: f64,
        y: f64,
        threshold_sq: f64,
        seed: Option<usize>,
    ) -> Option<(f64, usize)> {
        if self.segments.is_empty() {
            return None;
        }
        if self.segments.len() <= LINEAR_SCAN_MAX {
            let mut best = f64::INFINITY;
            let mut best_idx = 0usize;
            for seg_idx in 0..self.segments.len() {
                let d_sq = self.segments.point_segment_distance_sq(seg_idx, x, y);
                if d_sq <= threshold_sq {
                    return None;
                }
                if d_sq < best {
                    best = d_sq;
                    best_idx = seg_idx;
                }
            }
            return Some((best, best_idx));
        }
        self.bvh
            .nearest_distance_sq_above(&self.segments, x, y, threshold_sq, seed)
    }

    /// Indices of segments within `min_dist_sq + eps` of point (x, y).
    pub fn supporting_segments(&self, x: f64, y: f64, min_dist_sq: f64, eps: f64) -> Vec<usize> {
        let max_dist_sq = min_dist_sq + eps.abs().max(1e-14);
        if self.segments.is_empty() {
            return Vec::new();
        }
        if self.segments.len() <= LINEAR_SCAN_MAX {
            return linear_supporting_segments(&self.segments, x, y, max_dist_sq);
        }
        self.bvh
            .supporting_segments(&self.segments, x, y, max_dist_sq)
    }
}

fn linear_scan_nearest(segments: &SegmentIndex, x: f64, y: f64) -> (f64, usize) {
    let mut best_sq = f64::INFINITY;
    let mut best_idx = 0usize;
    for seg_idx in 0..segments.len() {
        let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
        if d_sq < best_sq {
            best_sq = d_sq;
            best_idx = seg_idx;
        }
    }
    (best_sq, best_idx)
}

fn linear_supporting_segments(
    segments: &SegmentIndex,
    x: f64,
    y: f64,
    max_dist_sq: f64,
) -> Vec<usize> {
    let mut supports = Vec::new();
    for seg_idx in 0..segments.len() {
        let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
        if d_sq <= max_dist_sq {
            supports.push(seg_idx);
        }
    }
    supports
}
