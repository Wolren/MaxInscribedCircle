use super::super::input::SegmentIndex;

/// Flat BVH (Bounding Volume Hierarchy) for nearest-boundary queries.
///
/// Provides accelerated nearest-distance and supporting-segment queries
/// over a set of indexed polygon segments.
#[derive(Debug, Clone)]
pub struct FlatBvh {
    segments_len: usize,
}

impl FlatBvh {
    /// Build a BVH from the given segment index.
    pub fn new(segments: &SegmentIndex) -> Self {
        Self {
            segments_len: segments.len(),
        }
    }

    /// Find the nearest segment to a point.
    /// Returns (squared_distance, segment_index).
    pub fn nearest_distance_sq(
        &self,
        segments: &SegmentIndex,
        x: f64,
        y: f64,
    ) -> Option<(f64, usize)> {
        if self.segments_len == 0 {
            return None;
        }
        let n = self.segments_len;
        let (d_sq, idx) = segments.batch_point_segment_distance_sq_with_index(x, y, 0, n);
        Some((d_sq, idx))
    }

    /// Find all segments within `max_dist_sq` of a point.
    pub fn supporting_segments(
        &self,
        segments: &SegmentIndex,
        x: f64,
        y: f64,
        max_dist_sq: f64,
    ) -> Vec<usize> {
        let mut supports = Vec::new();
        for seg_idx in 0..self.segments_len {
            let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
            if d_sq <= max_dist_sq {
                supports.push(seg_idx);
            }
        }
        supports
    }
}
