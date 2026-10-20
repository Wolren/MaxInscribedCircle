//! The BVH must answer exactly what a brute-force scan answers.
//!
//! Speed is only worth having if the index is indistinguishable from the scan
//! it replaces, so every case here compares the two directly.

use max_inscribed_circle::geo_types::{Coord, LineString, Polygon};
use max_inscribed_circle::mic::index::NearestBoundaryIndex;
use max_inscribed_circle::mic::input::{HostPolygon, SegmentIndex};

fn poly(pts: &[(f64, f64)]) -> Polygon<f64> {
    let mut v: Vec<Coord<f64>> = pts.iter().map(|p| Coord { x: p.0, y: p.1 }).collect();
    v.push(v[0]);
    Polygon::new(LineString::from(v), vec![])
}

fn regular(n: usize, r: f64) -> Polygon<f64> {
    let pts: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let a = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            (r * a.cos(), r * a.sin())
        })
        .collect();
    poly(&pts)
}

fn brute(seg: &SegmentIndex, x: f64, y: f64) -> (f64, usize) {
    let mut best = f64::INFINITY;
    let mut idx = 0;
    for i in 0..seg.len() {
        let d = seg.point_segment_distance_sq(i, x, y);
        if d < best {
            best = d;
            idx = i;
        }
    }
    (best, idx)
}

fn check(name: &str, p: &Polygon<f64>, probes: usize) {
    let host = HostPolygon::from_polygon(p).expect("host");
    let seg = SegmentIndex::from_host(&host);
    let nb = NearestBoundaryIndex::new(seg.clone());

    let (min_x, min_y, max_x, max_y) = host.bounds().expect("bounds");
    let mut checked = 0usize;
    for k in 0..probes {
        // Deterministic, spread over and beyond the box.
        let a = (k as f64) * 0.7548776662466927;
        let b = (k as f64) * 0.5698402909980532;
        let x = min_x - 0.1 * (max_x - min_x) + (max_x - min_x) * 1.2 * (a.fract());
        let y = min_y - 0.1 * (max_y - min_y) + (max_y - min_y) * 1.2 * (b.fract());

        let (want_d, _) = brute(&seg, x, y);
        let (got_d, got_i) = nb.nearest_distance_sq(x, y).expect("query");

        assert!(
            (want_d - got_d).abs() <= 1e-9 * want_d.max(1.0),
            "{name}: distance mismatch at ({x:.6},{y:.6}): brute {want_d:e} vs index {got_d:e}"
        );
        // Any segment at the same distance is a valid answer, so compare the
        // reported distance against the reported index rather than the index
        // against the scan's index.
        let via_idx = seg.point_segment_distance_sq(got_i, x, y);
        assert!(
            (via_idx - got_d).abs() <= 1e-9 * got_d.max(1.0),
            "{name}: index {got_i} does not carry the returned distance"
        );
        checked += 1;
    }
    assert!(checked == probes);
}

#[test]
fn bvh_matches_brute_force_on_many_shapes() {
    check(
        "square",
        &poly(&[(0.0, 0.0), (100.0, 0.0), (100.0, 100.0), (0.0, 100.0)]),
        400,
    );
    check(
        "triangle",
        &poly(&[(0.0, 0.0), (100.0, 0.0), (50.0, 86.0)]),
        400,
    );
    check(
        "L shape",
        &poly(&[
            (0.0, 0.0),
            (400.0, 0.0),
            (400.0, 1.0),
            (100.0, 1.0),
            (100.0, 300.0),
            (0.0, 300.0),
        ]),
        400,
    );
    check("regular 16-gon", &regular(16, 100.0), 400);
    check("regular 64-gon", &regular(64, 100.0), 400);
    check("regular 300-gon", &regular(300, 100.0), 600);

    // Polygon with a hole: the hole contributes real boundary segments.
    let frame = Polygon::new(
        LineString::from(vec![
            Coord { x: 0.0, y: 0.0 },
            Coord { x: 200.0, y: 0.0 },
            Coord { x: 200.0, y: 200.0 },
            Coord { x: 0.0, y: 200.0 },
            Coord { x: 0.0, y: 0.0 },
        ]),
        vec![LineString::from(vec![
            Coord { x: 60.0, y: 60.0 },
            Coord { x: 140.0, y: 60.0 },
            Coord { x: 140.0, y: 140.0 },
            Coord { x: 60.0, y: 140.0 },
            Coord { x: 60.0, y: 60.0 },
        ])],
    );
    check("frame with hole", &frame, 600);
}

#[test]
fn bvh_supporting_segments_matches_brute_force() {
    let p = regular(120, 50.0);
    let host = HostPolygon::from_polygon(&p).expect("host");
    let seg = SegmentIndex::from_host(&host);
    let nb = NearestBoundaryIndex::new(seg.clone());
    let (min_x, min_y, max_x, max_y) = host.bounds().expect("bounds");

    for k in 0..200 {
        let x = min_x + (max_x - min_x) * ((k as f64) * 0.61803398875).fract();
        let y = min_y + (max_y - min_y) * ((k as f64) * 0.41421356237).fract();
        let (d_sq, _) = nb.nearest_distance_sq(x, y).expect("query");
        let eps = d_sq.max(1.0) * 1e-10;
        let got = nb.supporting_segments(x, y, d_sq, eps);

        let max_sq = d_sq + eps;
        let want: Vec<usize> = (0..seg.len())
            .filter(|&i| seg.point_segment_distance_sq(i, x, y) <= max_sq)
            .collect();
        let mut got_sorted = got.clone();
        got_sorted.sort_unstable();
        let mut want_sorted = want.clone();
        want_sorted.sort_unstable();
        assert_eq!(
            got_sorted, want_sorted,
            "supporting set differs at ({x:.4},{y:.4})"
        );
    }
}

#[test]
fn bvh_handles_tiny_tables() {
    // Four segments: below the scan threshold, so this exercises the fallback
    // rather than the tree, and it must still answer exactly.
    let host =
        HostPolygon::from_polygon(&poly(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]))
            .expect("host");
    let nb = NearestBoundaryIndex::new(SegmentIndex::from_host(&host));
    let (d, _) = nb.nearest_distance_sq(5.0, 5.0).expect("query");
    assert!((d - 25.0).abs() < 1e-9, "expected 5^2, got {d}");
}

#[test]
fn bvh_builds_on_a_single_segment_table() {
    // A triangle has three non-degenerate segments, one short of the scan
    // threshold, so the tree is built but never consulted. Building it must
    // still be panic-free and the query must agree with the scan.
    let host =
        HostPolygon::from_polygon(&poly(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])).expect("host");
    let seg = SegmentIndex::from_host(&host);
    let nb = NearestBoundaryIndex::new(seg.clone());
    let (want, _) = brute(&seg, 0.2, 0.2);
    let (got, _) = nb.nearest_distance_sq(0.2, 0.2).expect("query");
    assert!((want - got).abs() <= 1e-12, "want {want:e} got {got:e}");
}
