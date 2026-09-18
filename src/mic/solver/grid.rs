//! Grid solver: tiles the bounding box with cells and refines the best one.
//!
//! The search is a branch and bound on one invariant: the distance to the
//! boundary is 1-Lipschitz, so for a cell of half-size `h` centred at `c`, no
//! point of that cell can be farther from the boundary than
//! `boundary_distance(c) + h * sqrt(2)`. That bound holds for every centre,
//! inside the polygon or not, which is what makes the pruning sound for cells
//! that straddle an edge. A cell leaves the queue only when its bound proves it
//! cannot beat the incumbent by more than the tolerance, so the queue drains
//! exactly when no region left in play can improve the answer. Nothing else
//! ends the search: there is no iteration budget, because a counter that stops
//! the loop early turns an unfinished search into an answer that looks
//! finished.
use std::collections::BinaryHeap;

use super::super::index::{NearestBoundaryIndex, PipIndex};
use super::super::input::HostPolygon;
use super::super::{MicApproach, MicResult, MicUsedEngine};
use crate::tuning;
use geo_types::Point;

const SQRT2: f64 = std::f64::consts::SQRT_2;

#[derive(Debug, Clone, Copy)]
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

/// Unsigned distance from a point to the nearest segment of the polygon.
///
/// This is the distance to the boundary as a set, so it is 1-Lipschitz
/// everywhere and bounds the clearance of every point of a cell around it,
/// whichever side of the polygon the queried point falls on.
#[inline]
fn boundary_distance(nb: &NearestBoundaryIndex, x: f64, y: f64) -> f64 {
    match nb.nearest_distance_sq(x, y) {
        Some((d, _)) if d.is_finite() && d > 0.0 => d.sqrt(),
        _ => 0.0,
    }
}

/// True when a cell of half-size `half` centred outside the polygon at `d`
/// from the boundary can hold no interior point.
///
/// If it held one, the segment from the inside point to the outside centre
/// would cross the boundary, and that crossing lies inside the cell, so the
/// centre would be at most `half * sqrt(2)` from the boundary. Finding the
/// centre farther away than that therefore proves the cell is empty of
/// interior points, and dropping it costs nothing.
#[inline]
fn empty_of_interior(inside: bool, d: f64, half: f64) -> bool {
    !inside && d > half * SQRT2
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
/// Tile the bounding box and refine the best cell into the inscribed circle.
///
/// The first level tiles the whole box, so every point of the polygon starts
/// inside some queued cell. Children wholly outside the box are dropped, which
/// is safe because a polygon is contained in its own bounding box. The radius
/// of the returned result is measured against the segment table at the returned
/// centre rather than read off the cell bound, so the reported circle always
/// fits inside the polygon.
pub fn solve_grid(
    host: &HostPolygon,
    pip: &PipIndex,
    nb: &NearestBoundaryIndex,
    tolerance: f64,
) -> Option<MicResult> {
    let (min_x, min_y, max_x, max_y) = host.bounds()?;
    let width = max_x - min_x;
    let height = max_y - min_y;
    let cell_size = width.max(height);

    if !(cell_size > 0.0) || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }

    // Half-size of a first-level cell: the box is covered by at most GRID_SIDE
    // cells along its longer axis and by however many are needed along the
    // shorter one, so a thin polygon is still covered end to end.
    let h_size = cell_size / (2.0 * tuning::MIC_GRID_SIDE as f64);
    let step = 2.0 * h_size;
    let cells_x = ((width / step).ceil() as usize).max(1);
    let cells_y = ((height / step).ceil() as usize).max(1);

    let mut queue: BinaryHeap<GridCell> = BinaryHeap::with_capacity(cells_x * cells_y * 2);
    let mut best_x = 0.0;
    let mut best_y = 0.0;
    let mut best_dist = 0.0;
    let mut found = false;

    for i in 0..cells_x {
        let cx = min_x + h_size + i as f64 * step;
        for j in 0..cells_y {
            let cy = min_y + h_size + j as f64 * step;
            let d = boundary_distance(nb, cx, cy);
            let inside = pip.contains_strict_xy(cx, cy);
            if empty_of_interior(inside, d, h_size) {
                continue;
            }
            if inside && d > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = d;
                found = true;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size,
                distance: if inside { d } else { 0.0 },
                max_dist: d + h_size * SQRT2,
            });
        }
    }

    if !found {
        return None;
    }

    let mut iterations = 0usize;

    while let Some(cell) = queue.pop() {
        iterations += 1;

        // The heap is ordered by the bound, so once the best remaining cell
        // cannot beat the incumbent by more than the tolerance, no cell can.
        let radius_bound = best_dist + tolerance;
        if cell.max_dist <= radius_bound {
            break;
        }

        let half = cell.h_size * 0.5;
        let children = [
            (cell.x - half, cell.y - half),
            (cell.x + half, cell.y - half),
            (cell.x - half, cell.y + half),
            (cell.x + half, cell.y + half),
        ];

        for (cx, cy) in children {
            // A cell that lies wholly outside the box holds no polygon point.
            if cx + half < min_x || cx - half > max_x || cy + half < min_y || cy - half > max_y {
                continue;
            }
            let d = boundary_distance(nb, cx, cy);
            let max_dist = d + half * SQRT2;
            // The bound test needs no containment result, so it runs first and
            // most cells never pay for one.
            if max_dist <= radius_bound {
                continue;
            }
            let inside = pip.contains_strict_xy(cx, cy);
            if empty_of_interior(inside, d, half) {
                continue;
            }
            if inside && d > best_dist {
                best_x = cx;
                best_y = cy;
                best_dist = d;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size: half,
                distance: if inside { d } else { 0.0 },
                max_dist,
            });
        }
    }

    // Measure, do not trust: the radius handed to the caller is the true
    // distance to the nearest segment of this polygon at this centre.
    let segments = nb.segments();
    let radius_sq = segments.batch_point_segment_distance_sq(best_x, best_y, 0, segments.len());
    if !radius_sq.is_finite() {
        return None;
    }
    let radius_sq = radius_sq.max(0.0);
    let radius = radius_sq.sqrt();

    let support_segments =
        nb.supporting_segments(best_x, best_y, radius_sq, radius_sq.max(1.0) * 1e-10);

    Some(MicResult {
        center: Point::new(best_x, best_y),
        radius,
        radius_sq,
        support_segments,
        candidate_count: iterations,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::GridSolver,
        component_index: None,
    })
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
/// True when the polygon has a single, strictly convex exterior ring.
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
/// True when any ring of the polygon is a hole.
pub fn has_holes(host: &HostPolygon) -> bool {
    host.rings.iter().any(|r| r.is_hole)
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
/// True when the polygon is convex and has no holes.
pub fn is_simple_shape(host: &HostPolygon) -> bool {
    is_convex_polygon(host) && !has_holes(host)
}
