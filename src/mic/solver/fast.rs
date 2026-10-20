//! Closed-form fast paths for the shape classes where the optimum is known.
//!
//! Each function recognises one class and returns the analytic circle for it,
//! or `None` when the ring is not that class. The pipeline re-measures every
//! result against the segment table before returning it, so a ring that merely
//! looks like the class falls through to the grid solver instead of producing a
//! circle it cannot support.

use super::super::input::HostPolygon;
use super::super::{MicApproach, MicResult, MicUsedEngine};
use crate::tuning::{MIC_BINARY_ITERS, MIC_EPS, MIC_EXPANSION_ITERS};
use geo_types::Point;
// ---------------------------------------------------------------------------
// Phase 0: Fast paths - analytical solvers, no workspace needed.
// ---------------------------------------------------------------------------

/// Triangle incenter - exact O(1) MIC. Trigger: ring_count==1, outer.len()==4.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn fast_triangle(host: &HostPolygon) -> Option<MicResult> {
    if host.ring_count() != 1 {
        return None;
    }
    let outer = host.outer_ring();
    if outer.len() != 4 {
        return None;
    }
    let a = outer[0];
    let b = outer[1];
    let c = outer[2];
    let la = (b[0] - c[0]).hypot(b[1] - c[1]);
    let lb = (a[0] - c[0]).hypot(a[1] - c[1]);
    let lc = (a[0] - b[0]).hypot(a[1] - b[1]);
    let perim = la + lb + lc;
    if perim <= 1e-14 {
        return None;
    }
    let cx = (la * a[0] + lb * b[0] + lc * c[0]) / perim;
    let cy = (la * a[1] + lb * b[1] + lc * c[1]) / perim;
    let area = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() * 0.5;
    let r = 2.0 * area / perim;
    if r <= 1e-6 {
        return None;
    }
    Some(MicResult {
        center: Point::new(cx, cy),
        radius: r,
        radius_sq: r * r,
        support_segments: vec![],
        candidate_count: 1,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::FastTriangle,
        component_index: None,
    })
}

/// Axis-aligned rectangle MIC - exact O(1). Trigger: ring_count==1, 4 vertices, right angles.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn fast_rectangle(host: &HostPolygon) -> Option<MicResult> {
    if host.ring_count() != 1 {
        return None;
    }
    let outer = host.outer_ring();
    if outer.len() != 5 {
        return None;
    }
    let v = [outer[0], outer[1], outer[2], outer[3]];
    // Check axis-aligned: all edges horizontal or vertical
    for i in 0..4 {
        let a = v[i];
        let b = v[(i + 1) % 4];
        let dx = (b[0] - a[0]).abs();
        let dy = (b[1] - a[1]).abs();
        if dx > 1e-12 && dy > 1e-12 {
            return None;
        }
        let len = dx + dy;
        if len <= 1e-12 {
            return None;
        }
    }
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in &v {
        min_x = min_x.min(p[0]);
        min_y = min_y.min(p[1]);
        max_x = max_x.max(p[0]);
        max_y = max_y.max(p[1]);
    }
    let width = max_x - min_x;
    let height = max_y - min_y;
    let r = width.min(height) * 0.5;
    if r <= 1e-12 {
        return None;
    }
    let cx = (min_x + max_x) * 0.5;
    let cy = (min_y + max_y) * 0.5;
    Some(MicResult {
        center: Point::new(cx, cy),
        radius: r,
        radius_sq: r * r,
        support_segments: vec![],
        candidate_count: 1,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::FastRectangle,
        component_index: None,
    })
}

/// Convex quadrilateral MIC via Chebyshev center in a local coordinate frame.
/// Preconditions: single ring, 4 vertices + closing point, convex and CCW.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn fast_convex_quad(host: &HostPolygon) -> Option<MicResult> {
    if host.ring_count() != 1 {
        return None;
    }
    let outer = host.outer_ring();
    if outer.len() != 5 {
        return None;
    }

    // Distinct vertices v0..v3 (ignore closing duplicate at index 4)
    let v = [outer[0], outer[1], outer[2], outer[3]];

    // Convexity check: all cross products same sign (CCW positive)
    let mut sign = 0i8;
    for i in 0..4 {
        let a = v[i];
        let b = v[(i + 1) % 4];
        let c = v[(i + 2) % 4];
        let cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if cross.abs() <= 1e-14 {
            return None;
        }
        let s = if cross > 0.0 { 1 } else { -1 };
        if sign == 0 {
            sign = s;
        } else if s != sign {
            return None;
        }
    }

    // --- Local frame: align edge v0->v1 with x-axis, y-axis points inward ---
    fn orthonormal_basis(p0: [f64; 2], p1: [f64; 2]) -> ((f64, f64), (f64, f64)) {
        let dx = p1[0] - p0[0];
        let dy = p1[1] - p0[1];
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 1e-14 {
            return ((1.0, 0.0), (0.0, 1.0));
        }
        let ux = dx / len;
        let uy = dy / len;
        // Inward normal for CCW polygon: rotate edge vector 90° clockwise -> (dy, -dx) but normalized as (-uy, ux)
        let vx = -uy;
        let vy = ux;
        ((ux, uy), (vx, vy))
    }
    fn to_local(p: [f64; 2], origin: [f64; 2], u: (f64, f64), v: (f64, f64)) -> (f64, f64) {
        let dx = p[0] - origin[0];
        let dy = p[1] - origin[1];
        (dx * u.0 + dy * u.1, dx * v.0 + dy * v.1)
    }
    fn from_local(x: f64, y: f64, origin: [f64; 2], u: (f64, f64), v: (f64, f64)) -> [f64; 2] {
        [origin[0] + x * u.0 + y * v.0, origin[1] + x * u.1 + y * v.1]
    }

    let origin = v[0];
    let (u, v_dir) = orthonormal_basis(v[0], v[1]);
    let q_local: [(f64, f64); 4] = [
        to_local(v[0], origin, u, v_dir),
        to_local(v[1], origin, u, v_dir),
        to_local(v[2], origin, u, v_dir),
        to_local(v[3], origin, u, v_dir),
    ];

    // Build edges: inward unit normals, line a*x + b*y + c = 0, signed distance = a*x+b*y+c
    #[derive(Clone, Copy)]
    struct Edge {
        a: f64,
        b: f64,
        c: f64,
    }
    let mut edges = [Edge {
        a: 0.0,
        b: 0.0,
        c: 0.0,
    }; 4];
    let orientation = sign as f64; // 1 for CCW, -1 for CW
    for i in 0..4 {
        let j = (i + 1) % 4;
        let (x0, y0) = q_local[i];
        let (x1, y1) = q_local[j];
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 1e-14 {
            return None;
        }
        // inward normal for CCW polygon: rotate (dx,dy) by +90° (left): (-dy, dx)
        // for CW polygon we negate the normal so it still points inward
        let a = -dy / len * orientation;
        let b = dx / len * orientation;
        // point on edge: a*x0 + b*y0 + c = 0  => c = -(a*x0 + b*y0)
        let c = -(a * x0 + b * y0);
        edges[i] = Edge { a, b, c };
    }

    // Upper bound on radius: half of polygon's smaller bbox dimension (safe overestimate)
    let bounds = host.bounds()?;
    let width = bounds.2 - bounds.0;
    let height = bounds.3 - bounds.1;
    let mut r_hi = width.min(height) * 0.5;
    let mut r_lo = 0.0;

    // Helper: min signed distance to any edge at (x,y)
    let min_dist = |eds: &[Edge; 4], x: f64, y: f64| -> f64 {
        eds.iter()
            .map(|e| e.a * x + e.b * y + e.c)
            .fold(f64::INFINITY, f64::min)
    };

    // Binary search for maximum feasible radius
    for _ in 0..MIC_BINARY_ITERS {
        let r_mid = 0.5 * (r_lo + r_hi);
        // Feasibility: exists point (x,y) such that for all edges: a*x+b*y+c >= r_mid
        let mut feasible = false;
        'search: for i in 0..4 {
            for j in (i + 1)..4 {
                let e1 = edges[i];
                let e2 = edges[j];
                // Solve: a1*x + b1*y = r_mid - c1,  a2*x + b2*y = r_mid - c2
                let a1 = e1.a;
                let b1 = e1.b;
                let rhs1 = r_mid - e1.c;
                let a2 = e2.a;
                let b2 = e2.b;
                let rhs2 = r_mid - e2.c;
                let det = a1 * b2 - a2 * b1;
                if det.abs() < MIC_EPS {
                    continue;
                }
                let inv = 1.0 / det;
                let x = (rhs1 * b2 - rhs2 * b1) * inv;
                let y = (a1 * rhs2 - a2 * rhs1) * inv;
                if !x.is_finite() || !y.is_finite() {
                    continue;
                }
                // Check all constraints
                for e in edges.iter() {
                    if e.a * x + e.b * y + e.c < r_mid - 1e-9 {
                        continue 'search;
                    }
                }
                feasible = true;
                break 'search;
            }
        }
        if feasible {
            r_lo = r_mid;
        } else {
            r_hi = r_mid;
        }
    }

    let r = r_lo;
    // Find a center point at radius r (intersection of two offset edges)
    let mut best_center = None;
    let mut best_min = -f64::INFINITY;
    for i in 0..4 {
        for j in (i + 1)..4 {
            let e1 = edges[i];
            let e2 = edges[j];
            let a1 = e1.a;
            let b1 = e1.b;
            let rhs1 = r - e1.c;
            let a2 = e2.a;
            let b2 = e2.b;
            let rhs2 = r - e2.c;
            let det = a1 * b2 - a2 * b1;
            if det.abs() < MIC_EPS {
                continue;
            }
            let inv = 1.0 / det;
            let x = (rhs1 * b2 - rhs2 * b1) * inv;
            let y = (a1 * rhs2 - a2 * rhs1) * inv;
            if !x.is_finite() || !y.is_finite() {
                continue;
            }
            let d_min = min_dist(&edges, x, y);
            if d_min > best_min {
                best_min = d_min;
                best_center = Some((x, y));
            }
        }
    }
    let (cx_local, cy_local) = best_center.unwrap_or_else(|| {
        // Fallback: centroid of vertices (should not happen)
        let cx = (q_local[0].0 + q_local[1].0 + q_local[2].0 + q_local[3].0) * 0.25;
        let cy = (q_local[0].1 + q_local[1].1 + q_local[2].1 + q_local[3].1) * 0.25;
        (cx, cy)
    });
    let center_world = from_local(cx_local, cy_local, origin, u, v_dir);

    Some(MicResult {
        center: Point::new(center_world[0], center_world[1]),
        radius: r,
        radius_sq: r * r,
        support_segments: vec![],
        candidate_count: 1,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::FastConvexQuad,
        component_index: None,
    })
}

/// General convex polygon MIC via Chebyshev center LP solved with binary search.
/// Works for any convex single-ring polygon with 3+ vertices.
/// For N edges, complexity is O(N^2 log precision) - much faster than candidate generation.
#[cfg_attr(feature = "hotpath", hotpath::measure)]
pub fn fast_convex_n(host: &HostPolygon) -> Option<MicResult> {
    if host.ring_count() != 1 {
        return None;
    }
    let outer = host.outer_ring();
    let n = outer.len().saturating_sub(1); // ignore closing point
    if n < 3 {
        return None;
    }

    // Convexity check + compute inward normals
    let mut sign = 0i8;
    let mut edges: Vec<(f64, f64, f64)> = Vec::with_capacity(n);
    for i in 0..n {
        let a = outer[i];
        let b = outer[(i + 1) % n];
        let c = outer[(i + 2) % n];
        let cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
        if cross.abs() > 1e-14 {
            let s = if cross > 0.0 { 1 } else { -1 };
            if sign == 0 {
                sign = s;
            } else if s != sign {
                return None;
            }
        }
        // inward normal for edge a->b
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let len = (dx * dx + dy * dy).sqrt();
        if len <= 1e-14 {
            return None;
        }
        let nx = -dy / len;
        let ny = dx / len;
        let nc = -(nx * a[0] + ny * a[1]);
        edges.push((nx, ny, nc));
    }

    // Flip normals if polygon is CW
    let orientation = sign as f64;
    if orientation < 0.0 {
        for e in edges.iter_mut() {
            e.0 = -e.0;
            e.1 = -e.1;
            e.2 = -e.2;
        }
    }

    let bounds = host.bounds()?;
    let width = bounds.2 - bounds.0;
    let height = bounds.3 - bounds.1;
    let mut r_hi = width.min(height) * 0.5;
    let mut r_lo = 0.0;

    // Binary search: feasibility = intersection of all half-planes a*x+b*y+c >= r is non-empty
    for _ in 0..MIC_EXPANSION_ITERS {
        let r_mid = 0.5 * (r_lo + r_hi);
        if convex_feasible(&edges, r_mid) {
            r_lo = r_mid;
        } else {
            r_hi = r_mid;
        }
    }

    let r = r_lo;
    // Extract a center point at radius r from best consecutive pair intersection
    let mut best_center = None;
    let mut best_min = -f64::INFINITY;
    for i in 0..n {
        let j = (i + 1) % n;
        let e1 = edges[i];
        let e2 = edges[j];
        let a1 = e1.0;
        let b1 = e1.1;
        let rhs1 = r - e1.2;
        let a2 = e2.0;
        let b2 = e2.1;
        let rhs2 = r - e2.2;
        let det = a1 * b2 - a2 * b1;
        if det.abs() < MIC_EPS {
            continue;
        }
        let inv = 1.0 / det;
        let x = (rhs1 * b2 - rhs2 * b1) * inv;
        let y = (a1 * rhs2 - a2 * rhs1) * inv;
        if !x.is_finite() || !y.is_finite() {
            continue;
        }
        let d_min = edges
            .iter()
            .map(|e| e.0 * x + e.1 * y + e.2)
            .fold(f64::INFINITY, f64::min);
        if d_min > best_min {
            best_min = d_min;
            best_center = Some((x, y));
        }
    }
    let (cx, cy) = best_center?;

    Some(MicResult {
        center: Point::new(cx, cy),
        radius: r,
        radius_sq: r * r,
        support_segments: vec![],
        candidate_count: 1,
        used_engine: MicUsedEngine::Grid,
        approach: MicApproach::FastConvexN,
        component_index: None,
    })
}

/// Check if intersection of all half-planes a*x+b*y+c >= r is non-empty for a convex polygon.
/// For convex polygons, it's sufficient to check consecutive edge-pair intersections.
fn convex_feasible(edges: &[(f64, f64, f64)], r: f64) -> bool {
    let n = edges.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let e1 = edges[i];
        let e2 = edges[j];
        let a1 = e1.0;
        let b1 = e1.1;
        let rhs1 = r - e1.2;
        let a2 = e2.0;
        let b2 = e2.1;
        let rhs2 = r - e2.2;
        let det = a1 * b2 - a2 * b1;
        if det.abs() < MIC_EPS {
            continue;
        }
        let inv = 1.0 / det;
        let x = (rhs1 * b2 - rhs2 * b1) * inv;
        let y = (a1 * rhs2 - a2 * rhs1) * inv;
        if !x.is_finite() || !y.is_finite() {
            continue;
        }
        let mut ok = true;
        for e in edges.iter() {
            if e.0 * x + e.1 * y + e.2 < r - 1e-9 {
                ok = false;
                break;
            }
        }
        if ok {
            return true;
        }
    }
    false
}
