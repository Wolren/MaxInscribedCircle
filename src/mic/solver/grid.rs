//! Grid solver: tiles the bounding box with cells and refines the best one.
//!
//! The search is a branch and bound on one invariant: the distance to the
//! boundary is 1-Lipschitz, so for a cell of half-size `h` centred at `c`, no
//! point of that cell can be farther from the boundary than
//! `boundary_distance(c) + h * sqrt(2)`. That bound holds for every centre,
//! inside the polygon or not, which is what makes the pruning sound for cells
//! that straddle an edge. A cell leaves the queue only when its bound proves it
//! cannot beat the incumbent by more than the tolerance, so the queue drains
//! exactly when no region left in play can improve the answer.
//!
//! The bound is equally loose in every direction, so a ridge of equal clearance
//! is resolved to `tolerance` along its whole length: cells grow like
//! `length / tolerance`. [`tuning::MIC_GRID_CELL_BUDGET`] caps that work, and an
//! exhausted budget stops the search and keeps the incumbent, which the hill
//! climb below then sharpens.
//!
//! The same looseness costs work on a broad plateau of near-optimal clearance.
//! A cell leaves the queue only once its bound proves it cannot reach the
//! incumbent, so a cell that could reach within a few percent must be refined a
//! long way first. On the vendored 290-parcel corpus that is the dominant cost,
//! and it follows from the Lipschitz bound rather than from this loop.
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
    max_dist: f64,
    /// Segment that achieves the nearest distance at this cell's centre, used to
    /// seed the next level's distance query so it starts from a near hit.
    support: usize,
}

impl PartialEq for GridCell {
    fn eq(&self, other: &Self) -> bool {
        self.max_dist == other.max_dist
    }
}

impl Eq for GridCell {}

impl PartialOrd for GridCell {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
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

    if cell_size <= 0.0 || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }

    // Half-size of a first-level cell: the box is covered by at most GRID_SIDE
    // cells along its longer axis and by however many are needed along the
    // shorter one, so a thin polygon is still covered end to end.
    //
    // The seed density is capped by how expensive one boundary query is. A
    // 25x25 seed is 625 queries, which is more than the search itself pops on
    // every shape measured, and on a 256-gon those queries cost 394 us of a
    // 597 us solve. A coarse seed costs nothing in accuracy because the bound
    // still prunes and the incumbent is measured at the end; it only costs a
    // few extra refinement levels, each of which is also just queries.
    let side = tuning::MIC_GRID_SIDE;
    let h_size = cell_size / (2.0 * side as f64);
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
                max_dist: d + h_size * SQRT2,
                support: usize::MAX,
            });
        }
    }

    let mut iterations = 0usize;
    // True when the loop stopped on the budget instead of on the bound, which
    // means the incumbent may still be short of the peak.
    let mut exhausted_budget = false;

    while let Some(cell) = queue.pop() {
        iterations += 1;

        // A ridge of equal clearance costs length / tolerance cells, so the
        // budget stops shapes the bisection cannot finish in reasonable time.
        // Stopping keeps the incumbent: every centre the search accepts is
        // inside the polygon and its radius is measured against the segment
        // table at the end, so the circle returned below still fits.
        if iterations > tuning::MIC_GRID_CELL_BUDGET {
            exhausted_budget = true;
            break;
        }

        // The heap is ordered by the bound, so once the best remaining cell
        // cannot beat the incumbent by more than the prune gap, no cell can.
        // With no incumbent yet there is nothing to compare against, so the
        // search continues and the first interior centre becomes the incumbent.
        let radius_bound = best_dist + tolerance;
        if found && cell.max_dist <= radius_bound {
            break;
        }

        let half = cell.h_size * 0.5;
        let children = [
            (cell.x - half, cell.y - half),
            (cell.x + half, cell.y - half),
            (cell.x - half, cell.y + half),
            (cell.x + half, cell.y + half),
        ];

        // A child is pruned when `d + half*sqrt(2) <= best + tolerance`, which
        // rearranges to `d <= best - max(0, half*sqrt(2) - tolerance)`. That
        // right-hand side is the threshold below which the exact distance no
        // longer matters: any child centre closer than it is discarded either
        // way. Asking the index for "the nearest distance, or a value under
        // this bar" lets the descent stop at the first segment inside the bar
        // instead of proving the true minimum. This is the technique polylabel
        // uses, and it makes the common hopeless child far cheaper.
        let slack = half * SQRT2 - (radius_bound - best_dist);
        let threshold = if found && slack > 0.0 {
            Some(best_dist - slack)
        } else {
            None
        };

        for (cx, cy) in children {
            // A cell that lies wholly outside the box holds no polygon point.
            if cx + half < min_x || cx - half > max_x || cy + half < min_y || cy - half > max_y {
                continue;
            }
            let (d, support) = match threshold {
                Some(t) => match nb.nearest_distance_sq_above(
                    cx,
                    cy,
                    (t.max(0.0)) * (t.max(0.0)),
                    Some(cell.support),
                ) {
                    Some(v) => v,
                    // Proven within the bar: this child cannot improve.
                    None => continue,
                },
                None => {
                    let Some(v) = nb.nearest_distance_sq(cx, cy) else {
                        continue;
                    };
                    (v.0, v.1)
                }
            };
            let d = d.max(0.0).sqrt();
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
                // The seed row can miss the polygon entirely when the shape is
                // thinner than the first cell. The refinement is what finds the
                // interior then, so an interior centre seen here counts as an
                // incumbent: without this the bound test never engages and the
                // search refines the whole box down to the tolerance.
                found = true;
            }
            queue.push(GridCell {
                x: cx,
                y: cy,
                h_size: half,
                max_dist,
                support,
            });
        }
    }

    // Measure, do not trust: the radius handed to the caller is the true
    // distance to the nearest segment of this polygon at this centre. Without
    // an incumbent there is no centre to measure, and the seeds covered the
    // box, so the only way to get here is a polygon with no interior point at
    // all, which the caller reports as an error.
    if !found {
        return None;
    }

    // Polish only when the search ran out of budget. A search that stopped on
    // the bound is already within `tolerance` of the optimum, so climbing it
    // again spends probes to recover nothing. A ridge, by contrast, outruns the
    // budget before the bisection closes and can sit a percent short, which is
    // the case this walk exists for. It only ever accepts a centre that is
    // strictly better and strictly inside.
    let (best_x, best_y) = if exhausted_budget {
        let (px, py, _) = polish(nb, pip, best_x, best_y, best_dist, tolerance, h_size);
        (px, py)
    } else {
        (best_x, best_y)
    };

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

/// Climb the clearance field from an incumbent centre until no probe improves.
///
/// Each round probes the eight compass directions at a shrinking radius, and
/// accepts a probe only when it is strictly inside the polygon and strictly
/// better. The radius halves whenever a round finds nothing, so the walk stops
/// once the step falls under `tol`. Returns the best centre found and its
/// clearance, never worse than the input.
fn polish(
    nb: &NearestBoundaryIndex,
    pip: &PipIndex,
    start_x: f64,
    start_y: f64,
    start_d: f64,
    tol: f64,
    initial_step: f64,
) -> (f64, f64, f64) {
    let mut best_x = start_x;
    let mut best_y = start_y;
    let mut best_d = start_d;
    let mut step = initial_step;

    // A ridge peak sits where two supports balance, so a compass walk reaches
    // it in a handful of rounds; the cap only guards against a pathological
    // field that keeps improving by a vanishing amount.
    for _ in 0..256 {
        if step <= tol || !step.is_finite() {
            break;
        }
        let mut moved = false;
        for (dx, dy) in [
            (-1.0, 0.0),
            (1.0, 0.0),
            (0.0, -1.0),
            (0.0, 1.0),
            (-1.0, -1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (1.0, 1.0),
        ] {
            let nx = best_x + dx * step;
            let ny = best_y + dy * step;
            let Some((d_sq, _)) = nb.nearest_distance_sq(nx, ny) else {
                continue;
            };
            let d = d_sq.max(0.0).sqrt();
            if d > best_d && pip.contains_strict_xy(nx, ny) {
                best_x = nx;
                best_y = ny;
                best_d = d;
                moved = true;
            }
        }
        if !moved {
            step *= 0.5;
        }
    }
    (best_x, best_y, best_d)
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
