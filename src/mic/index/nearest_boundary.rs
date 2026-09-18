//! Nearest-boundary distance queries over a segment table, backed by a
//! uniform grid or, when requested, a flat BVH.
use super::super::input::SegmentIndex;
use super::bvh::FlatBvh;

const GRID_MIN_CELLS: usize = 4;
const GRID_MAX_CELLS: usize = 128;

#[derive(Debug, Clone)]
struct GridIndex {
    cells: Vec<Vec<usize>>,
    nx: usize,
    ny: usize,
    bbox_min_x: f64,
    bbox_min_y: f64,
    cell_w: f64,
    cell_h: f64,
}

impl GridIndex {
    fn from_segments(segments: &SegmentIndex) -> Self {
        if segments.is_empty() {
            return Self {
                cells: Vec::new(),
                nx: 0,
                ny: 0,
                bbox_min_x: 0.0,
                bbox_min_y: 0.0,
                cell_w: 1.0,
                cell_h: 1.0,
            };
        }

        // Compute overall bbox (sequential - cheap)
        let mut bbox_min_x = f64::INFINITY;
        let mut bbox_min_y = f64::INFINITY;
        let mut bbox_max_x = f64::NEG_INFINITY;
        let mut bbox_max_y = f64::NEG_INFINITY;

        for idx in 0..segments.len() {
            bbox_min_x = bbox_min_x.min(segments.bbox_minx[idx]);
            bbox_min_y = bbox_min_y.min(segments.bbox_miny[idx]);
            bbox_max_x = bbox_max_x.max(segments.bbox_maxx[idx]);
            bbox_max_y = bbox_max_y.max(segments.bbox_maxy[idx]);
        }

        let span_x = (bbox_max_x - bbox_min_x).max(1e-12);
        let span_y = (bbox_max_y - bbox_min_y).max(1e-12);

        let target_cells = (segments.len() as f64 * 0.5).sqrt().ceil() as usize;
        let nx = target_cells.clamp(GRID_MIN_CELLS, GRID_MAX_CELLS);
        let ny = target_cells.clamp(GRID_MIN_CELLS, GRID_MAX_CELLS);
        let cell_w = span_x / nx as f64;
        let cell_h = span_y / ny as f64;

        // Build grid cells sequentially
        let mut cells: Vec<Vec<usize>> = vec![Vec::new(); nx * ny];
        for seg_idx in 0..segments.len() {
            let min_ci = ((segments.bbox_minx[seg_idx] - bbox_min_x) / cell_w).floor() as isize;
            let max_ci = ((segments.bbox_maxx[seg_idx] - bbox_min_x) / cell_w).ceil() as isize;
            let min_cj = ((segments.bbox_miny[seg_idx] - bbox_min_y) / cell_h).floor() as isize;
            let max_cj = ((segments.bbox_maxy[seg_idx] - bbox_min_y) / cell_h).ceil() as isize;

            let min_ci = min_ci.max(0).min((nx - 1) as isize);
            let max_ci = max_ci.max(0).min((nx - 1) as isize);
            let min_cj = min_cj.max(0).min((ny - 1) as isize);
            let max_cj = max_cj.max(0).min((ny - 1) as isize);

            for cj in min_cj..=max_cj {
                for ci in min_ci..=max_ci {
                    let cell_idx = cj as usize * nx + ci as usize;
                    cells[cell_idx].push(seg_idx);
                }
            }
        }

        Self {
            cells,
            nx,
            ny,
            bbox_min_x,
            bbox_min_y,
            cell_w,
            cell_h,
        }
    }

    fn cell_of(&self, x: f64, y: f64) -> Option<(isize, isize)> {
        if x < self.bbox_min_x || y < self.bbox_min_y {
            return None;
        }
        let ci = ((x - self.bbox_min_x) / self.cell_w) as isize;
        let cj = ((y - self.bbox_min_y) / self.cell_h) as isize;
        if ci < 0 || cj < 0 || ci >= self.nx as isize || cj >= self.ny as isize {
            return None;
        }
        Some((ci, cj))
    }

    fn cell_bbox_dist_sq(&self, x: f64, y: f64, ci: isize, cj: isize) -> f64 {
        let cell_min_x = self.bbox_min_x + ci as f64 * self.cell_w;
        let cell_max_x = cell_min_x + self.cell_w;
        let cell_min_y = self.bbox_min_y + cj as f64 * self.cell_h;
        let cell_max_y = cell_min_y + self.cell_h;

        let dx = if x < cell_min_x {
            cell_min_x - x
        } else if x > cell_max_x {
            x - cell_max_x
        } else {
            0.0
        };
        let dy = if y < cell_min_y {
            cell_min_y - y
        } else if y > cell_max_y {
            y - cell_max_y
        } else {
            0.0
        };
        dx * dx + dy * dy
    }

    #[inline]
    fn nearest_via_grid(&self, segments: &SegmentIndex, x: f64, y: f64) -> Option<(f64, usize)> {
        let (cx, cy) = self.cell_of(x, y)?;

        let mut best_sq = f64::INFINITY;
        let mut best_idx = 0usize;

        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let max_ring = nx.max(ny);

        for ring in 0..=max_ring {
            let mut ring_min_dist = f64::INFINITY;

            let ci0 = (cx - ring).max(0);
            let ci1 = (cx + ring).min(nx - 1);
            let cj0 = (cy - ring).max(0);
            let cj1 = (cy + ring).min(ny - 1);

            if ci0 > ci1 || cj0 > cj1 {
                break;
            }

            if ring == 0 {
                let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx, cy);
                ring_min_dist = cell_dist_sq;
                let cell_idx = (cy * nx + cx) as usize;
                for &seg_idx in &self.cells[cell_idx] {
                    let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                    if d_sq < best_sq {
                        best_sq = d_sq;
                        best_idx = seg_idx;
                    }
                }
            } else {
                if cj0 <= cy - ring && cy - ring <= cj1 {
                    for ci in ci0..=ci1 {
                        if (ci - cx).abs().max(ring) == ring {
                            let cell_dist_sq = self.cell_bbox_dist_sq(x, y, ci, cy - ring);
                            ring_min_dist = ring_min_dist.min(cell_dist_sq);
                            if cell_dist_sq <= best_sq {
                                let cell_idx = ((cy - ring) * nx + ci) as usize;
                                for &seg_idx in &self.cells[cell_idx] {
                                    let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                    if d_sq < best_sq {
                                        best_sq = d_sq;
                                        best_idx = seg_idx;
                                    }
                                }
                            }
                        }
                    }
                }
                if cj0 <= cy + ring && cy + ring <= cj1 {
                    for ci in ci0..=ci1 {
                        if (ci - cx).abs().max(ring) == ring {
                            let cell_dist_sq = self.cell_bbox_dist_sq(x, y, ci, cy + ring);
                            ring_min_dist = ring_min_dist.min(cell_dist_sq);
                            if cell_dist_sq <= best_sq {
                                let cell_idx = ((cy + ring) * nx + ci) as usize;
                                for &seg_idx in &self.cells[cell_idx] {
                                    let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                    if d_sq < best_sq {
                                        best_sq = d_sq;
                                        best_idx = seg_idx;
                                    }
                                }
                            }
                        }
                    }
                }
                if ci0 <= cx - ring && cx - ring <= ci1 {
                    for cj in cj0..=cj1 {
                        if ring.max((cj - cy).abs()) == ring {
                            let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx - ring, cj);
                            ring_min_dist = ring_min_dist.min(cell_dist_sq);
                            if cell_dist_sq <= best_sq {
                                let cell_idx = (cj * nx + cx - ring) as usize;
                                for &seg_idx in &self.cells[cell_idx] {
                                    let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                    if d_sq < best_sq {
                                        best_sq = d_sq;
                                        best_idx = seg_idx;
                                    }
                                }
                            }
                        }
                    }
                }
                if ci0 <= cx + ring && cx + ring <= ci1 {
                    for cj in cj0..=cj1 {
                        if ring.max((cj - cy).abs()) == ring {
                            let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx + ring, cj);
                            ring_min_dist = ring_min_dist.min(cell_dist_sq);
                            if cell_dist_sq <= best_sq {
                                let cell_idx = (cj * nx + cx + ring) as usize;
                                for &seg_idx in &self.cells[cell_idx] {
                                    let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                    if d_sq < best_sq {
                                        best_sq = d_sq;
                                        best_idx = seg_idx;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if ring_min_dist >= best_sq {
                break;
            }
        }

        Some((best_sq, best_idx))
    }

    fn supporting_via_grid(
        &self,
        segments: &SegmentIndex,
        x: f64,
        y: f64,
        max_dist_sq: f64,
    ) -> Vec<usize> {
        let Some((cx, cy)) = self.cell_of(x, y) else {
            return linear_supporting_segments(segments, x, y, max_dist_sq);
        };

        let mut supports = Vec::new();
        let nx = self.nx as isize;
        let ny = self.ny as isize;
        let max_ring = nx.max(ny);

        for ring in 0..=max_ring {
            let mut ring_min_dist = f64::INFINITY;

            let ci0 = (cx - ring).max(0);
            let ci1 = (cx + ring).min(nx - 1);
            let cj0 = (cy - ring).max(0);
            let cj1 = (cy + ring).min(ny - 1);

            if ci0 > ci1 || cj0 > cj1 {
                break;
            }

            if ring == 0 {
                let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx, cy);
                ring_min_dist = cell_dist_sq;
                if cell_dist_sq <= max_dist_sq {
                    let cell_idx = (cy * nx + cx) as usize;
                    for &seg_idx in &self.cells[cell_idx] {
                        let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                        if d_sq <= max_dist_sq {
                            supports.push(seg_idx);
                        }
                    }
                }
            } else {
                if cj0 <= cy - ring && cy - ring <= cj1 {
                    for ci in ci0..=ci1 {
                        let cell_dist_sq = self.cell_bbox_dist_sq(x, y, ci, cy - ring);
                        ring_min_dist = ring_min_dist.min(cell_dist_sq);
                        if cell_dist_sq <= max_dist_sq {
                            let cell_idx = ((cy - ring) * nx + ci) as usize;
                            for &seg_idx in &self.cells[cell_idx] {
                                let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                if d_sq <= max_dist_sq {
                                    supports.push(seg_idx);
                                }
                            }
                        }
                    }
                }
                if cj0 <= cy + ring && cy + ring <= cj1 {
                    for ci in ci0..=ci1 {
                        let cell_dist_sq = self.cell_bbox_dist_sq(x, y, ci, cy + ring);
                        ring_min_dist = ring_min_dist.min(cell_dist_sq);
                        if cell_dist_sq <= max_dist_sq {
                            let cell_idx = ((cy + ring) * nx + ci) as usize;
                            for &seg_idx in &self.cells[cell_idx] {
                                let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                if d_sq <= max_dist_sq {
                                    supports.push(seg_idx);
                                }
                            }
                        }
                    }
                }
                if ci0 <= cx - ring && cx - ring <= ci1 {
                    for cj in cj0..=cj1 {
                        let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx - ring, cj);
                        ring_min_dist = ring_min_dist.min(cell_dist_sq);
                        if cell_dist_sq <= max_dist_sq {
                            let cell_idx = (cj * nx + cx - ring) as usize;
                            for &seg_idx in &self.cells[cell_idx] {
                                let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                if d_sq <= max_dist_sq {
                                    supports.push(seg_idx);
                                }
                            }
                        }
                    }
                }
                if ci0 <= cx + ring && cx + ring <= ci1 {
                    for cj in cj0..=cj1 {
                        let cell_dist_sq = self.cell_bbox_dist_sq(x, y, cx + ring, cj);
                        ring_min_dist = ring_min_dist.min(cell_dist_sq);
                        if cell_dist_sq <= max_dist_sq {
                            let cell_idx = (cj * nx + cx + ring) as usize;
                            for &seg_idx in &self.cells[cell_idx] {
                                let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
                                if d_sq <= max_dist_sq {
                                    supports.push(seg_idx);
                                }
                            }
                        }
                    }
                }
            }

            if ring_min_dist > max_dist_sq {
                break;
            }
        }

        supports
    }
}

#[inline]
fn linear_scan_nearest(
    segments: &SegmentIndex,
    x: f64,
    y: f64,
    start_seg: usize,
) -> Option<(f64, usize)> {
    let n = segments.len();
    if n >= 8 {
        // Use SIMD batch for larger segment counts
        let (d_sq, idx) = segments.batch_point_segment_distance_sq_with_index(x, y, start_seg, n);
        return Some((d_sq, idx));
    }

    let mut best_sq = f64::INFINITY;
    let mut best_idx = start_seg;

    for seg_idx in 0..n {
        let bbox_lb = point_to_bbox_distance_sq(
            x,
            y,
            segments.bbox_minx[seg_idx],
            segments.bbox_miny[seg_idx],
            segments.bbox_maxx[seg_idx],
            segments.bbox_maxy[seg_idx],
        );
        if bbox_lb > best_sq {
            continue;
        }
        let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
        if d_sq < best_sq {
            best_sq = d_sq;
            best_idx = seg_idx;
        }
    }

    Some((best_sq, best_idx))
}

fn linear_supporting_segments(
    segments: &SegmentIndex,
    x: f64,
    y: f64,
    max_dist_sq: f64,
) -> Vec<usize> {
    let mut supports = Vec::new();
    for seg_idx in 0..segments.len() {
        let bbox_lb = point_to_bbox_distance_sq(
            x,
            y,
            segments.bbox_minx[seg_idx],
            segments.bbox_miny[seg_idx],
            segments.bbox_maxx[seg_idx],
            segments.bbox_maxy[seg_idx],
        );
        if bbox_lb > max_dist_sq {
            continue;
        }
        let d_sq = segments.point_segment_distance_sq(seg_idx, x, y);
        if d_sq <= max_dist_sq {
            supports.push(seg_idx);
        }
    }
    supports
}

/// Nearest-boundary distance queries over segment table, accelerated with a uniform grid.
#[derive(Debug, Clone)]
pub struct NearestBoundaryIndex {
    segments: SegmentIndex,
    grid: Option<GridIndex>,
    bvh: Option<FlatBvh>,
}

impl NearestBoundaryIndex {
    /// Build the index over `segments`, using a BVH when `use_bvh` is true.
    pub fn new(segments: SegmentIndex, use_bvh: bool) -> Self {
        if use_bvh {
            let bvh = FlatBvh::new(&segments);
            Self {
                segments,
                grid: None,
                bvh: Some(bvh),
            }
        } else {
            let grid = GridIndex::from_segments(&segments);
            Self {
                segments,
                grid: Some(grid),
                bvh: None,
            }
        }
    }

    #[inline]
    /// The segment table the index was built from.
    pub fn segments(&self) -> &SegmentIndex {
        &self.segments
    }

    #[inline]
    /// Squared distance to the closest boundary and its segment index, if any.
    pub fn nearest_distance_sq(&self, x: f64, y: f64) -> Option<(f64, usize)> {
        if let Some(bvh) = &self.bvh {
            return bvh.nearest_distance_sq(&self.segments, x, y);
        }
        let grid = self.grid.as_ref().unwrap();
        if self.segments.is_empty() {
            return None;
        }
        let n = self.segments.len();
        // Use SIMD batch when segment count is reasonable (faster than grid for small polygons)
        #[cfg(feature = "simd")]
        if n <= 256 {
            return Some(
                self.segments
                    .batch_point_segment_distance_sq_with_index(x, y, 0, n),
            );
        }
        #[cfg(not(feature = "simd"))]
        if n <= 32 {
            return linear_scan_nearest(&self.segments, x, y, 0);
        }
        if grid.nx == 0 {
            return linear_scan_nearest(&self.segments, x, y, 0);
        }
        grid.nearest_via_grid(&self.segments, x, y)
    }

    /// Indices of segments within `min_dist_sq + eps` of point (x, y).
    pub fn supporting_segments(&self, x: f64, y: f64, min_dist_sq: f64, eps: f64) -> Vec<usize> {
        let max_dist_sq = min_dist_sq + eps.abs().max(1e-14);
        if let Some(bvh) = &self.bvh {
            return bvh.supporting_segments(&self.segments, x, y, max_dist_sq);
        }
        let grid = self.grid.as_ref().unwrap();
        if grid.nx == 0 {
            return linear_supporting_segments(&self.segments, x, y, max_dist_sq);
        }
        grid.supporting_via_grid(&self.segments, x, y, max_dist_sq)
    }
}

fn point_to_bbox_distance_sq(
    x: f64,
    y: f64,
    min_x: f64,
    min_y: f64,
    max_x: f64,
    max_y: f64,
) -> f64 {
    let dx = if x < min_x {
        min_x - x
    } else if x > max_x {
        x - max_x
    } else {
        0.0
    };
    let dy = if y < min_y {
        min_y - y
    } else if y > max_y {
        y - max_y
    } else {
        0.0
    };
    dx * dx + dy * dy
}
