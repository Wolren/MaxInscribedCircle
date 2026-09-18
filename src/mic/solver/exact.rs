//! Exact MIC solver: candidate enumeration over boundary vertices,
//! segment triples, ear circumcenters, and Voronoi sites.
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use geo::Centroid;
use geo_types::Point;
use rustc_hash::FxHashSet;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use super::super::index::NearestBoundaryIndex;
use super::super::input::{HostPolygon, SegmentIndex};
use super::super::workspace::{MicCandidate, MicWorkspace};
use super::super::{MicApproach, MicError, MicOptions, MicResult, MicUsedEngine, RobustMode};
use crate::tuning::{
    MIC_BASE_SEGS_PER_RING, MIC_BASE_SS_SEG_CAP, MIC_BASE_SS_VERT_CAP, MIC_BASE_TRIPLE_CAP,
    MIC_BINARY_ITERS, MIC_CANDIDATE_QUANTIZE, MIC_EPS, MIC_EXPANSION_ITERS, MIC_EXT_SEGS_PER_RING,
    MIC_EXT_SS_SEG_CAP, MIC_EXT_SS_VERT_CAP, MIC_EXT_TRIPLE_CAP,
};

const CANDIDATE_QUANTIZE: f64 = MIC_CANDIDATE_QUANTIZE;

// Baseline caps for simple (convex, few segments) polygons.
const BASE_TRIPLE_CAP: usize = MIC_BASE_TRIPLE_CAP;
const BASE_SS_SEG_CAP: usize = MIC_BASE_SS_SEG_CAP;
const BASE_SS_VERT_CAP: usize = MIC_BASE_SS_VERT_CAP;
const BASE_SEGS_PER_RING: usize = MIC_BASE_SEGS_PER_RING;

const EXT_TRIPLE_CAP: usize = MIC_EXT_TRIPLE_CAP;
const EXT_SS_SEG_CAP: usize = MIC_EXT_SS_SEG_CAP;
const EXT_SS_VERT_CAP: usize = MIC_EXT_SS_VERT_CAP;
const EXT_SEGS_PER_RING: usize = MIC_EXT_SEGS_PER_RING;

/// Check if a vertex at `vert_idx` in ring `ring_idx` is reflex (interior angle > 180°).
fn reflex_vertex_in_ring(host: &HostPolygon, ring_idx: usize, vert_idx: usize) -> bool {
    let meta = &host.rings[ring_idx];
    let coords = &host.coords[meta.start..meta.end];
    let n = if coords.len() >= 2 && coords.first() == coords.last() {
        coords.len() - 1
    } else {
        coords.len()
    };
    if n < 3 {
        return false;
    }
    let idx = vert_idx % n;
    let prev = coords[(idx + n - 1) % n];
    let cur = coords[idx];
    let next = coords[(idx + 1) % n];
    let cross = (cur[0] - prev[0]) * (next[1] - cur[1]) - (cur[1] - prev[1]) * (next[0] - cur[0]);
    if meta.is_hole {
        cross > 1e-14
    } else {
        cross < -1e-14
    }
}

/// Count reflex vertices across all rings (concavity measure).
fn count_reflex_vertices(host: &HostPolygon, seg_index: &SegmentIndex) -> usize {
    use std::collections::HashSet;
    let mut reflex_set = HashSet::new();
    for seg_idx in 0..seg_index.len() {
        let rid = seg_index.ring_id[seg_idx];
        let eid = seg_index.edge_id[seg_idx];
        if reflex_vertex_in_ring(host, rid, eid) || reflex_vertex_in_ring(host, rid, eid + 1) {
            reflex_set.insert((rid, eid));
        }
    }
    reflex_set.len()
}

/// Compute adaptive cap tier based on polygon complexity.
/// Complex polygons (reflex vertices, holes, high segment count) get extended caps
/// to ensure the candidate set covers the true MIC's support segments.
fn caps_for(
    seg_count: usize,
    hole_count: usize,
    reflex_count: usize,
) -> (usize, usize, usize, usize) {
    // Trigger extended caps if any of:
    //   - >3 holes (cross-ring MIC constraints)
    //   - >8 reflex vertices (deep concavity)
    //   - >200 segments (dense boundary)
    let complex = hole_count > 3 || reflex_count > 8 || seg_count > 200;
    if complex {
        (
            EXT_TRIPLE_CAP,
            EXT_SS_SEG_CAP,
            EXT_SS_VERT_CAP,
            EXT_SEGS_PER_RING,
        )
    } else {
        (
            BASE_TRIPLE_CAP,
            BASE_SS_SEG_CAP,
            BASE_SS_VERT_CAP,
            BASE_SEGS_PER_RING,
        )
    }
}

fn quantize_origin(host: &HostPolygon) -> (f64, f64) {
    let Some((min_x, min_y, max_x, _max_y)) = host.bounds() else {
        return (0.0, 0.0);
    };
    let span_x = (max_x - min_x).max(1.0);
    (min_x - span_x * 0.1, min_y - span_x * 0.1)
}

#[cfg_attr(feature = "hotpath", hotpath::measure)]
/// Run the exact solver over a prepared workspace and return the best circle.
pub fn solve_exact(
    workspace: &mut MicWorkspace,
    opts: &MicOptions,
) -> std::result::Result<MicResult, MicError> {
    workspace.clear_candidates();
    let mut seen = FxHashSet::<(i64, i64)>::default();
    let q_origin = quantize_origin(&workspace.host);

    // 1. Centroid
    if let Some(c) = workspace.host.polygon.centroid() {
        push_candidate(
            &mut workspace.candidate_buf,
            &mut seen,
            c.x(),
            c.y(),
            q_origin,
        );
    }

    // 2. All boundary vertices
    let vertices = workspace.host.unique_vertices();
    for v in &vertices {
        push_candidate(
            &mut workspace.candidate_buf,
            &mut seen,
            v[0],
            v[1],
            q_origin,
        );
    }

    // 3. Segment midpoints
    for seg_idx in 0..workspace.nb_index.segments().len() {
        let (mx, my) = workspace.nb_index.segments().midpoint(seg_idx);
        push_candidate(&mut workspace.candidate_buf, &mut seen, mx, my, q_origin);
    }

    // Compute adaptive caps based on polygon complexity.
    let hole_count = workspace.host.rings.iter().filter(|r| r.is_hole).count();
    let reflex_count = count_reflex_vertices(&workspace.host, &workspace.nb_index.segments());
    let seg_count = workspace.nb_index.segments().len();
    let (triple_cap, ss_seg_cap, ss_vert_cap, segs_per_ring) =
        caps_for(seg_count, hole_count, reflex_count);

    // 4. Segment-triple incenters - reflex-biased sampling (Gap B)
    generate_segment_triple_candidates(
        &workspace.nb_index.segments(),
        &workspace.host,
        &mut seen,
        &mut workspace.candidate_buf,
        q_origin,
        triple_cap,
        segs_per_ring,
    );

    // 5. CDT circumcenters - only for polygons with enough vertices
    // where cheaper generators may miss the optimum.
    let total_vertices = workspace.host.coords.len();
    if total_vertices > 30 {
        generate_cdt_candidates(
            &workspace.host,
            &mut seen,
            &mut workspace.candidate_buf,
            q_origin,
        );
    }

    // 6. Ear circumcenters - ALL rings including holes
    generate_ear_candidates_all_rings(
        &workspace.host,
        &mut seen,
        &mut workspace.candidate_buf,
        q_origin,
    );

    // 7. Filtered: seg-seg-vertex bisector candidates + vertex-triple circumcenters
    if matches!(opts.robust_mode, RobustMode::Filtered) {
        let lines = precompute_segment_lines(&workspace.nb_index.segments());
        generate_ssv_candidates(
            &workspace.nb_index.segments(),
            &lines,
            &vertices,
            &mut seen,
            &mut workspace.candidate_buf,
            q_origin,
            ss_seg_cap,
            ss_vert_cap,
        );

        let sampled = sample_vertices(&vertices, 48);
        for i in 0..sampled.len() {
            for j in i + 1..sampled.len() {
                for k in j + 1..sampled.len() {
                    if let Some((cx, cy)) = circumcenter(sampled[i], sampled[j], sampled[k]) {
                        push_candidate(&mut workspace.candidate_buf, &mut seen, cx, cy, q_origin);
                    }
                }
            }
        }
    }

    if workspace.candidate_buf.is_empty() {
        return Err(MicError::NoCircleFound);
    }

    let candidate_count = workspace.candidate_buf.len();
    let pip_index = &workspace.pip_index;
    let nb_index = &workspace.nb_index;
    let mut best_idx: Option<usize> = None;
    let mut best_radius_sq = 0.0;

    let candidate_buf = &mut workspace.candidate_buf;

    for (i, cand) in candidate_buf.iter_mut().enumerate() {
        if !pip_index.contains_strict_xy(cand.x, cand.y) {
            continue;
        }
        let Some((radius_sq, _nearest_idx)) = nb_index.nearest_distance_sq(cand.x, cand.y) else {
            continue;
        };
        if !radius_sq.is_finite() || radius_sq <= 0.0 {
            continue;
        }
        cand.radius_sq = radius_sq;

        if cand.radius_sq > best_radius_sq {
            best_radius_sq = cand.radius_sq;
            best_idx = Some(i);
        }
    }

    let best = best_idx
        .map(|i| candidate_buf[i].clone())
        .ok_or(MicError::NoCircleFound)?;

    // Phase 4: Quadtree refinement from best candidate.
    // Closes the gap between discrete candidate optimum and true MIC.
    let (min_x, min_y, max_x, max_y) = workspace.host.bounds().unwrap_or((0.0, 0.0, 1.0, 1.0));
    let diameter = (max_x - min_x).hypot(max_y - min_y).max(1.0);
    let refine_tol = diameter * 1e-6;
    let seed_h = best.radius_sq.sqrt().max(1e-12) * 0.5;
    let (ref_x, ref_y, ref_r) = quadtree_search(
        best.x,
        best.y,
        best.radius_sq.sqrt(),
        seed_h,
        pip_index,
        nb_index,
        refine_tol,
        100,
    );

    let center = Point::new(ref_x, ref_y);
    let ref_r_sq = ref_r * ref_r;
    let support_eps = ref_r_sq.max(1.0) * 1e-10;
    let support_segments = nb_index.supporting_segments(ref_x, ref_y, ref_r_sq, support_eps);

    Ok(MicResult {
        center,
        radius: ref_r,
        radius_sq: ref_r_sq,
        support_segments,
        candidate_count,
        used_engine: MicUsedEngine::Exact,
        approach: MicApproach::VoronoiExact,
        component_index: None,
    })
}

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
        used_engine: MicUsedEngine::Exact,
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
        used_engine: MicUsedEngine::Exact,
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
        used_engine: MicUsedEngine::Exact,
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
        used_engine: MicUsedEngine::Exact,
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

// ---------------------------------------------------------------------------
// Quantization helpers
// ---------------------------------------------------------------------------

#[inline]
fn push_candidate(
    buf: &mut Vec<MicCandidate>,
    seen: &mut FxHashSet<(i64, i64)>,
    x: f64,
    y: f64,
    q_origin: (f64, f64),
) {
    if !x.is_finite() || !y.is_finite() {
        return;
    }
    let qx = quantize(x - q_origin.0);
    let qy = quantize(y - q_origin.1);
    if !seen.insert((qx, qy)) {
        return;
    }
    buf.push(MicCandidate {
        x,
        y,
        radius_sq: 0.0,
    });
}

#[inline]
fn quantize(v: f64) -> i64 {
    (v * CANDIDATE_QUANTIZE).round() as i64
}

fn sample_vertices(vertices: &[[f64; 2]], max_vertices: usize) -> Vec<[f64; 2]> {
    if vertices.len() <= max_vertices {
        return vertices.to_vec();
    }
    let step = ((vertices.len() as f64) / (max_vertices as f64)).ceil() as usize;
    vertices.iter().step_by(step.max(1)).copied().collect()
}

/// Sample segments with ring awareness - guarantees at least MIN_SEGS_PER_RING
/// from each ring before distributing the remaining budget by segment count.
fn sample_segments_ring_aware(
    seg_index: &SegmentIndex,
    max_total: usize,
    min_per_ring: usize,
) -> Vec<usize> {
    let n = seg_index.len();
    if n <= max_total || n == 0 {
        return (0..n).collect();
    }
    if seg_index.ring_id.is_empty() {
        return Vec::new();
    }

    // Find ring boundaries in the flat segment list
    let mut ring_starts: Vec<usize> = vec![0];
    let mut last_rid = seg_index.ring_id[0];
    for i in 1..n {
        if seg_index.ring_id[i] != last_rid {
            ring_starts.push(i);
            last_rid = seg_index.ring_id[i];
        }
    }
    let num_rings = ring_starts.len();
    let mut ring_ends = ring_starts[1..].to_vec();
    ring_ends.push(n);

    // Allocate: min_per_ring guaranteed, remainder by segment count
    let guaranteed = min_per_ring * num_rings;
    let remaining = if max_total > guaranteed {
        max_total - guaranteed
    } else {
        0
    };

    let mut result = Vec::with_capacity(max_total);
    for ri in 0..num_rings {
        let start = ring_starts[ri];
        let end = ring_ends[ri];
        let count = end - start;
        let alloc = if count <= min_per_ring {
            count
        } else {
            let extra = remaining * count / n;
            min_per_ring + extra.min(count - min_per_ring)
        };
        if count <= alloc {
            for idx in start..end {
                result.push(idx);
            }
        } else {
            let step = count / alloc;
            for i in (0..count).step_by(step.max(1)).take(alloc) {
                result.push(start + i);
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Segment-triple incenter candidates (segment-segment-segment Voronoi vertices)
// ---------------------------------------------------------------------------

fn generate_segment_triple_candidates(
    seg_index: &SegmentIndex,
    host: &HostPolygon,
    seen: &mut FxHashSet<(i64, i64)>,
    candidate_buf: &mut Vec<MicCandidate>,
    q_origin: (f64, f64),
    triple_cap: usize,
    _segs_per_ring: usize,
) {
    let n = seg_index.len();
    if n < 3 {
        return;
    }
    let lines = precompute_segment_lines(seg_index);

    // Reflex-biased sampling (Gap B): unconditionally include segments
    // adjacent to reflex vertices, then fill remaining budget with
    // ring-aware uniform sampling from non-reflex segments.
    let sampled = if n <= triple_cap {
        (0..n).collect()
    } else {
        let mut result: Vec<usize> = (0..n)
            .filter(|&si| {
                reflex_vertex_in_ring(host, seg_index.ring_id[si], seg_index.edge_id[si])
                    || reflex_vertex_in_ring(host, seg_index.ring_id[si], seg_index.edge_id[si] + 1)
            })
            .collect();
        result.sort_unstable();
        result.dedup();
        if result.len() >= triple_cap {
            result.truncate(triple_cap)
        } else {
            let reflex_set: FxHashSet<usize> = result.iter().copied().collect();
            let remaining = triple_cap - result.len();
            let non_reflex: Vec<usize> = (0..n).filter(|i| !reflex_set.contains(i)).collect();
            if !non_reflex.is_empty() {
                let step = (non_reflex.len() / remaining).max(1);
                result.extend(non_reflex.iter().step_by(step).take(remaining).copied());
            }
        }
        result
    };
    for ii in 0..sampled.len() {
        let i = sampled[ii];
        for jj in ii + 1..sampled.len() {
            let j = sampled[jj];
            for kk in jj + 1..sampled.len() {
                let k = sampled[kk];
                if let Some((cx, cy)) = segment_incenter(&lines, i, j, k) {
                    push_candidate(candidate_buf, seen, cx, cy, q_origin);
                }
            }
        }
    }
}

fn segment_incenter(lines: &[SegmentLine], i: usize, j: usize, k: usize) -> Option<(f64, f64)> {
    let li = &lines[i];
    let lj = &lines[j];
    let lk = &lines[k];
    let a_x = li.nx - lj.nx;
    let a_y = li.ny - lj.ny;
    let d_ij = li.c - lj.c;
    let b_x = li.nx - lk.nx;
    let b_y = li.ny - lk.ny;
    let d_ik = li.c - lk.c;
    let det = a_x * b_y - a_y * b_x;
    if det.abs() <= 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    let x = (d_ij * b_y - d_ik * a_y) * inv_det;
    let y = (a_x * d_ik - b_x * d_ij) * inv_det;
    if !x.is_finite() || !y.is_finite() {
        return None;
    }
    let d_i = li.nx * x + li.ny * y - li.c;
    if d_i <= 0.0 {
        return None;
    }
    Some((x, y))
}

struct SegmentLine {
    nx: f64,
    ny: f64,
    c: f64,
}

fn precompute_segment_lines(seg_index: &SegmentIndex) -> Vec<SegmentLine> {
    let mut lines = Vec::with_capacity(seg_index.len());
    for idx in 0..seg_index.len() {
        let dx = seg_index.dir_x[idx];
        let dy = seg_index.dir_y[idx];
        let len = seg_index.len_sq[idx].sqrt().max(1e-300);
        let inv_len = 1.0 / len;
        let is_hole = seg_index.is_hole_edge[idx];
        // Outer ring (CCW): inward = rotate left = (-dy, dx)
        // Hole ring (CW): inward = rotate right = (dy, -dx)
        let (nx, ny) = if !is_hole {
            (-dy * inv_len, dx * inv_len)
        } else {
            (dy * inv_len, -dx * inv_len)
        };
        let c = nx * seg_index.ax[idx] + ny * seg_index.ay[idx];
        lines.push(SegmentLine { nx, ny, c });
    }
    lines
}

// ---------------------------------------------------------------------------
// Seg-seg-vertex bisector candidates (Filtered mode, small budget)
// Catches the rare segment-segment-vertex Voronoi vertices not covered by
// segment-triple (3 segments) or CDT (3 vertices).
// ---------------------------------------------------------------------------

fn generate_ssv_candidates(
    seg_index: &SegmentIndex,
    lines: &[SegmentLine],
    vertices: &[[f64; 2]],
    seen: &mut FxHashSet<(i64, i64)>,
    candidate_buf: &mut Vec<MicCandidate>,
    q_origin: (f64, f64),
    ss_seg_cap: usize,
    ss_vert_cap: usize,
) {
    if seg_index.len() < 2 || vertices.is_empty() {
        return;
    }

    let sampled_segs = sample_segments_ring_aware(seg_index, ss_seg_cap, 2);
    let max_verts = ss_vert_cap.min(vertices.len());
    let sampled_verts: Vec<[f64; 2]> = if vertices.len() <= max_verts {
        vertices.to_vec()
    } else {
        let step = vertices.len() / max_verts;
        vertices
            .iter()
            .step_by(step.max(1))
            .copied()
            .take(max_verts)
            .collect()
    };

    for ii in 0..sampled_segs.len() {
        let i = sampled_segs[ii];
        let li = &lines[i];
        for jj in ii + 1..sampled_segs.len() {
            let j = sampled_segs[jj];
            let lj = &lines[j];

            let nx = li.nx - lj.nx;
            let ny = li.ny - lj.ny;
            let n_len_sq = nx * nx + ny * ny;
            if n_len_sq <= 1e-14 {
                continue;
            }
            let d_ij = li.c - lj.c;
            let inv_n2 = 1.0 / n_len_sq;
            let c0x = nx * d_ij * inv_n2;
            let c0y = ny * d_ij * inv_n2;
            let n_len = n_len_sq.sqrt();
            let dx = -ny / n_len;
            let dy = nx / n_len;
            let dist0 = li.nx * c0x + li.ny * c0y - li.c;
            let nd = li.nx * dx + li.ny * dy;
            let coeff_a = 1.0 - nd * nd;
            if coeff_a.abs() <= 1e-14 {
                continue;
            }
            let inv_2a = 0.5 / coeff_a;

            for v in &sampled_verts {
                let dvx = c0x - v[0];
                let dvy = c0y - v[1];
                let delta_sq = dvx * dvx + dvy * dvy;
                let delta_dot_d = dvx * dx + dvy * dy;
                let coeff_b = 2.0 * (delta_dot_d - dist0 * nd);
                let coeff_c = delta_sq - dist0 * dist0;
                let disc = coeff_b * coeff_b - 4.0 * coeff_a * coeff_c;
                if disc < 0.0 {
                    continue;
                }
                let sqrt_disc = disc.sqrt();
                for t in [
                    (-coeff_b + sqrt_disc) * inv_2a,
                    (-coeff_b - sqrt_disc) * inv_2a,
                ] {
                    let cx = c0x + t * dx;
                    let cy = c0y + t * dy;
                    if !cx.is_finite() || !cy.is_finite() {
                        continue;
                    }
                    if li.nx * cx + li.ny * cy - li.c <= 0.0 {
                        continue;
                    }
                    if lj.nx * cx + lj.ny * cy - lj.c <= 0.0 {
                        continue;
                    }
                    push_candidate(candidate_buf, seen, cx, cy, q_origin);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// CDT Voronoi vertex generator (constrained Delaunay triangulation)
// ---------------------------------------------------------------------------

fn generate_cdt_candidates(
    host: &HostPolygon,
    seen: &mut FxHashSet<(i64, i64)>,
    candidate_buf: &mut Vec<MicCandidate>,
    q_origin: (f64, f64),
) {
    let mut cdt: ConstrainedDelaunayTriangulation<Point2<f64>> =
        ConstrainedDelaunayTriangulation::new();

    for ring in &host.rings {
        let coords = &host.coords[ring.start..ring.end];
        let n = if coords.len() >= 2 && coords.first() == coords.last() {
            coords.len() - 1
        } else {
            coords.len()
        };
        if n < 3 {
            continue;
        }

        // Insert vertices
        let mut handles = Vec::with_capacity(n);
        for i in 0..n {
            if let Ok(h) = cdt.insert(Point2::new(coords[i][0], coords[i][1])) {
                handles.push(h);
            }
        }
        // Add constraint edges. Crossing edges are split at the intersection
        // instead of panicking: the identity vertex constructor keeps the
        // split vertices at the exact intersection positions.
        for i in 0..handles.len() {
            let j = (i + 1) % handles.len();
            let _ = cdt.add_constraint_and_split(handles[i], handles[j], |p| p);
        }
    }

    for face in cdt.inner_faces() {
        let verts = face.vertices();
        let a = [verts[0].position().x, verts[0].position().y];
        let b = [verts[1].position().x, verts[1].position().y];
        let c = [verts[2].position().x, verts[2].position().y];
        if let Some((cx, cy)) = circumcenter(a, b, c) {
            push_candidate(candidate_buf, seen, cx, cy, q_origin);
        }
    }
}

// ---------------------------------------------------------------------------
// Ear circumcenter candidates - ALL rings including holes
// ---------------------------------------------------------------------------

fn generate_ear_candidates_all_rings(
    host: &HostPolygon,
    seen: &mut FxHashSet<(i64, i64)>,
    candidate_buf: &mut Vec<MicCandidate>,
    q_origin: (f64, f64),
) {
    for ring in &host.rings {
        let coords = &host.coords[ring.start..ring.end];
        let n = if coords.len() >= 2 && coords.first() == coords.last() {
            coords.len() - 1
        } else {
            coords.len()
        };
        if n < 3 {
            continue;
        }

        let verts = &coords[..n];
        let is_hole = ring.is_hole;

        for i in 0..n {
            let prev = if i == 0 { n - 1 } else { i - 1 };
            let next = if i + 1 >= n { 0 } else { i + 1 };
            let a = verts[prev];
            let b = verts[i];
            let c = verts[next];
            let cross = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
            // Outer ring CCW: convex when cross > 0
            // Hole ring CW: convex when cross < 0
            let is_convex = if !is_hole {
                cross > 1e-14
            } else {
                cross < -1e-14
            };
            if is_convex {
                if let Some((cx, cy)) = circumcenter(a, b, c) {
                    push_candidate(candidate_buf, seen, cx, cy, q_origin);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 4: Quadtree refinement pass
// ---------------------------------------------------------------------------

const SQRT_2: f64 = std::f64::consts::SQRT_2;

struct QuadCell(f64, f64, f64, f64);

impl PartialEq for QuadCell {
    fn eq(&self, other: &Self) -> bool {
        self.3 == other.3
    }
}
impl Eq for QuadCell {}
impl PartialOrd for QuadCell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.3.partial_cmp(&other.3).map(|o| o.reverse())
    }
}
impl Ord for QuadCell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.3
            .partial_cmp(&other.3)
            .unwrap_or(Ordering::Equal)
            .reverse()
    }
}

/// Run quadtree refinement with configurable iteration budget.
/// `seed_h`: initial cell half-side. For refinement (tight seed, ~best_r/2).
/// For primary solve (wide seed, ~diameter/4).
fn quadtree_search(
    init_x: f64,
    init_y: f64,
    init_r: f64,
    seed_h: f64,
    pip: &super::super::index::PipIndex,
    nb: &NearestBoundaryIndex,
    tol: f64,
    max_iters: usize,
) -> (f64, f64, f64) {
    let mut queue: BinaryHeap<QuadCell> = BinaryHeap::new();
    queue.push(QuadCell(init_x, init_y, seed_h, init_r + seed_h * SQRT_2));
    let mut best_x = init_x;
    let mut best_y = init_y;
    let mut best_r = init_r;
    let mut iters = 0usize;
    while let Some(QuadCell(cx, cy, h, upper)) = queue.pop() {
        iters += 1;
        if iters > max_iters {
            break;
        }
        if upper <= best_r + tol {
            break;
        }
        if h < tol {
            continue;
        }
        let h2 = h * 0.5;
        for (dx, dy) in [(-h2, -h2), (h2, -h2), (-h2, h2), (h2, h2)] {
            let nx = cx + dx;
            let ny = cy + dy;
            if !pip.contains_strict_xy(nx, ny) {
                continue;
            }
            let Some((r2, _)) = nb.nearest_distance_sq(nx, ny) else {
                continue;
            };
            let r = r2.sqrt();
            if r > best_r {
                best_r = r;
                best_x = nx;
                best_y = ny;
            }
            let ub = r + h2 * SQRT_2;
            if ub > best_r + tol {
                queue.push(QuadCell(nx, ny, h2, ub));
            }
        }
    }
    (best_x, best_y, best_r)
}

// ---------------------------------------------------------------------------
// Circumcenter of three 2D points
// ---------------------------------------------------------------------------

fn circumcenter(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> Option<(f64, f64)> {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    if d.abs() <= 1e-14 {
        return None;
    }
    let a2 = a[0] * a[0] + a[1] * a[1];
    let b2 = b[0] * b[0] + b[1] * b[1];
    let c2 = c[0] * c[0] + c[1] * c[1];
    let ux = (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d;
    let uy = (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d;
    if ux.is_finite() && uy.is_finite() {
        Some((ux, uy))
    } else {
        None
    }
}
