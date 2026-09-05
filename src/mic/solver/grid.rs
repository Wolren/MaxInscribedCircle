use std::collections::BinaryHeap;

use super::super::index::{NearestBoundaryIndex, PipIndex};
use super::super::input::HostPolygon;
use super::super::{MicApproach, MicResult, MicUsedEngine};
use geo_types::Point;

const SQRT2: f64 = std::f64::consts::SQRT_2;

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct GridCell {
    x: f64,
    y: f64,
    h_size: f64,
    distance: f64,
    max_dist: f64,
}

impl PartialEq for GridCell {
    fn eq(&self, other: &Self) -> bool {
        self.max_dist == other.max_dist
    }
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.max_dist.partial_cmp(&other.max_dist)
    }
}

impl Ord for GridCell {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.max_dist
            .partial_cmp(&other.max_dist)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn solve_grid(
    host: &HostPolygon,
    pip: &PipIndex,
    nb: &NearestBoundaryIndex,
    tolerance: f64,
) -> Option<MicResult> {
    let bounds = host.bounds()?;
    let (min_x, min_y, max_x, max_y) = bounds;

    let width = max_x - min_x;
    let height = max_y - min_y;
    let cell_size = width.max(height);

    if cell_size == 0.0 {
        return None;
    }

    let grid_side = 25;
    let h_size = cell_size / (grid_side as f64);

    // Pre-allocate priority queue with estimated capacity to reduce reallocations
    let initial_cell_count = grid_side * grid_side;
    let mut queue: BinaryHeap<GridCell> = BinaryHeap::with_capacity(initial_cell_count * 2);

    // Precompute base offsets to reduce arithmetic in inner loop
    let base_x = min_x + 0.5 * h_size;
    let base_y = min_y + 0.5 * h_size;

    for i in 0..grid_side {
        let cx_base = base_x + i as f64 * h_size;
        for j in 0..grid_side {
            let cx = cx_base;
            let cy = base_y + j as f64 * h_size;

            // Fast overall bbox reject: point outside polygon's AABB cannot be inside
            if cx > max_x || cy > max_y {
                continue;
            }

            // PIP check - skips most cells
            if !pip.contains_strict_xy(cx, cy) {
                continue;
            }

            // Nearest boundary distance
            let Some((dist_sq, _)) = nb.nearest_distance_sq(cx, cy) else {
                continue;
            };
            if !dist_sq.is_finite() || dist_sq <= 0.0 {
                continue;
            }

            let distance = dist_sq.sqrt();
            let max_dist = distance + h_size * SQRT2;

            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size,
                distance,
                max_dist,
            });
        }
    }

    if queue.is_empty() {
        return None;
    }

    let mut farthest_x = 0.0;
    let mut farthest_y = 0.0;
    let mut farthest_dist = 0.0;
    let mut farthest_dist_sq = 0.0;
    let mut found_initial = false;

    let max_iterations = compute_max_iterations(host, tolerance);
    let mut iterations = 0;

    while let Some(cell) = queue.pop() {
        iterations += 1;
        if iterations > max_iterations {
            break;
        }

        if found_initial && cell.max_dist <= farthest_dist {
            break;
        }

        if !pip.contains_strict_xy(cell.x, cell.y) {
            continue;
        }

        let Some((dist_sq, _)) = nb.nearest_distance_sq(cell.x, cell.y) else {
            continue;
        };

        if !dist_sq.is_finite() || dist_sq <= 0.0 {
            continue;
        }

        let dist = dist_sq.sqrt();

        // Update farthest cell if this one is better
        if !found_initial || dist > farthest_dist {
            farthest_x = cell.x;
            farthest_y = cell.y;
            farthest_dist = dist;
            farthest_dist_sq = dist_sq;
            found_initial = true;
        }

        if cell.h_size > tolerance {
            let half_h = cell.h_size * 0.5;
            // Split into 4 sub-cells (GEOS-style)
            let sub_cells = [
                (cell.x - half_h, cell.y - half_h),
                (cell.x + half_h, cell.y - half_h),
                (cell.x - half_h, cell.y + half_h),
                (cell.x + half_h, cell.y + half_h),
            ];

            for (nx, ny) in sub_cells {
                // Quick overall bbox reject - avoids PIP for clearly exterior points
                if nx < min_x || nx > max_x || ny < min_y || ny > max_y {
                    continue;
                }
                // PIP check
                if !pip.contains_strict_xy(nx, ny) {
                    continue;
                }
                // Nearest boundary distance
                let Some((dist_sq, _)) = nb.nearest_distance_sq(nx, ny) else {
                    continue;
                };
                if !dist_sq.is_finite() || dist_sq <= 0.0 {
                    continue;
                }
                let distance = dist_sq.sqrt();
                let max_dist = distance + half_h * SQRT2;
                // Only push if it could potentially improve current best
                if max_dist > farthest_dist + tolerance {
                    queue.push(GridCell {
                        x: nx,
                        y: ny,
                        h_size: half_h,
                        distance,
                        max_dist,
                    });
                }
            }
        }
    }

    if !found_initial {
        return None;
    }

    let support_segments = nb.supporting_segments(
        farthest_x,
        farthest_y,
        farthest_dist_sq,
        farthest_dist_sq.max(1.0) * 1e-10,
    );

    Some(MicResult {
        center: Point::new(farthest_x, farthest_y),
        radius: farthest_dist,
        radius_sq: farthest_dist_sq,
        support_segments,
        candidate_count: iterations,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::GridSolver,
        component_index: None,
    })
}
fn compute_max_iterations(host: &HostPolygon, tolerance: f64) -> usize {
    let Some(bounds) = host.bounds() else {
        return 10000;
    };
    let (min_x, min_y, max_x, max_y) = bounds;
    let width = max_x - min_x;
    let height = max_y - min_y;
    let area = width * height;

    let tolerance_sq = tolerance * tolerance;
    if tolerance_sq <= 0.0 {
        return 10000;
    }

    let cell_area = tolerance_sq;
    let cell_count = (area / cell_area).ceil() as usize;

    let safety_factor = 10;
    (cell_count * safety_factor).min(100000).max(100)
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn is_convex_polygon(host: &HostPolygon) -> bool {
    if host.ring_count() > 1 {
        return false;
    }

    let outer = host.outer_ring();
    if outer.len() < 4 {
        return false;
    }

    let n = outer.len() - 1;
    let mut sign = 0i8;

    for i in 0..n {
        let p = outer[(i + n - 1) % n];
        let c = outer[i];
        let nxt = outer[(i + 1) % n];

        let cross = (c[0] - p[0]) * (nxt[1] - c[1]) - (c[1] - p[1]) * (nxt[0] - c[0]);
        if cross.abs() <= 1e-14 {
            return false;
        }

        let s = if cross > 0.0 { 1 } else { -1 };
        if sign == 0 {
            sign = s;
        } else if s != sign {
            return false;
        }
    }

    true
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn has_holes(host: &HostPolygon) -> bool {
    host.rings.iter().any(|r| r.is_hole)
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn is_simple_shape(host: &HostPolygon) -> bool {
    is_convex_polygon(host) && !has_holes(host)
}
