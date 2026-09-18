//! C ABI for the MIC solvers (holes + multipolygon).
//!
//! Every exported function returns 0 on success, -1 on error, and never
//! lets a Rust panic cross into the host process. Ported from IGE `ige-c`
//! with the `ige_` prefix replaced by `mic_`. The GEOS fallback path stays
//! compiled out unless the `geos` feature is enabled.

use super::mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon, MicEngine, MicOptions,
    MicUsedEngine, RobustMode,
};
use geo_types::{Coord, LineString, Polygon};
use libc::{c_double, c_int, size_t};
use std::slice;

/// One hole ring: flat coords + length.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct HoleInput {
    /// Flat [x0, y0, x1, y1, ...] hole ring coordinates, or null when unused.
    pub coords: *const c_double,
    /// Number of doubles in `coords`; even and at least 6.
    pub coords_len: size_t,
}

/// MIC result.
#[repr(C)]
pub struct MicResultC {
    /// Circle center x coordinate.
    pub center_x: c_double,
    /// Circle center y coordinate.
    pub center_y: c_double,
    /// Circle radius in input coordinate units.
    pub radius: c_double,
    /// Squared radius, radius * radius.
    pub radius_sq: c_double,
    /// 0 = exact, 1 = geos fallback, 2 = grid.
    pub used_engine: c_int,
    /// Number of candidate centers evaluated during the solve.
    pub candidate_count: size_t,
    /// -1 when not a multipolygon solve.
    pub component_index: c_int,
}

/// MIC options.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct MicOptionsC {
    /// 0 = ExactOnly, 1 = FallbackOnly, 2 = ExactThenGeos.
    pub engine: c_int,
    /// 0 = FastF64, 1 = Filtered.
    pub robust_mode: c_int,
}

impl Default for MicOptionsC {
    fn default() -> Self {
        Self {
            engine: 2,
            robust_mode: 1,
        }
    }
}

/// One multipolygon component: exterior + holes.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct MicPolyInput {
    /// Flat [x0, y0, x1, y1, ...] exterior ring coordinates.
    pub coords: *const c_double,
    /// Number of doubles in `coords`; even and at least 6.
    pub coords_len: size_t,
    /// Pointer to `holes_len` hole rings, or null for no holes.
    pub holes: *const HoleInput,
    /// Number of entries in `holes`; ignored when `holes` is null.
    pub holes_len: size_t,
}

/// Decode a flat [x0, y0, x1, y1, ...] array into coords. Validates length.
///
/// # Safety
/// `coords` must point to `coords_len` readable doubles, or be null
/// (null returns None).
pub unsafe fn decode_coords(
    coords: *const c_double,
    coords_len: size_t,
) -> Option<Vec<Coord<f64>>> {
    if coords.is_null() {
        return None;
    }
    if coords_len < 6 || !coords_len.is_multiple_of(2) {
        return None;
    }
    let s = slice::from_raw_parts(coords, coords_len);
    let mut out = Vec::with_capacity(coords_len / 2);
    for i in (0..coords_len).step_by(2) {
        let (x, y) = (s[i], s[i + 1]);
        if !x.is_finite() || !y.is_finite() {
            return None;
        }
        out.push(Coord { x, y });
    }
    Some(out)
}

/// Decode exterior + optional holes into a polygon.
///
/// # Safety
/// Same contract as [`decode_coords`]; `holes` must point to `holes_len`
/// readable entries, or be null (null means no holes).
pub unsafe fn decode_polygon(
    coords: *const c_double,
    coords_len: size_t,
    holes: *const HoleInput,
    holes_len: size_t,
) -> Option<Polygon<f64>> {
    let ext = decode_coords(coords, coords_len)?;
    let mut interiors = Vec::new();
    if !holes.is_null() {
        let hs = slice::from_raw_parts(holes, holes_len);
        for h in hs {
            let ring = decode_coords(h.coords, h.coords_len)?;
            if ring.len() < 3 {
                return None;
            }
            interiors.push(LineString::from(ring));
        }
    }
    Some(Polygon::new(LineString::from(ext), interiors))
}

/// Run a closure behind the panic boundary. Panic becomes -1.
pub fn guard(f: impl FnOnce() -> c_int + std::panic::UnwindSafe) -> c_int {
    std::panic::catch_unwind(f).unwrap_or(-1)
}

fn engine_from(v: c_int) -> MicEngine {
    match v {
        0 => MicEngine::ExactOnly,
        1 => MicEngine::FallbackOnly,
        _ => MicEngine::ExactThenGeos,
    }
}

fn robust_from(v: c_int) -> RobustMode {
    match v {
        0 => RobustMode::FastF64,
        _ => RobustMode::Filtered,
    }
}

fn engine_to(v: MicUsedEngine) -> c_int {
    match v {
        MicUsedEngine::Exact => 0,
        MicUsedEngine::GeosFallback => 1,
        MicUsedEngine::Grid => 2,
    }
}

fn fill(r: &super::mic::MicResult, out: *mut MicResultC) -> c_int {
    if out.is_null() {
        return -1;
    }
    unsafe {
        *out = MicResultC {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            used_engine: engine_to(r.used_engine),
            candidate_count: r.candidate_count as size_t,
            component_index: r.component_index.map(|x| x as c_int).unwrap_or(-1),
        };
    }
    0
}

fn opts_from(raw: MicOptionsC) -> MicOptions {
    MicOptions {
        engine: engine_from(raw.engine),
        robust_mode: robust_from(raw.robust_mode),
        use_bvh: false,
        use_parallel: false,
    }
}

/// Solve MIC. Holes optional (null = none).
///
/// Returns 0 on success, -1 on error. `FallbackOnly` without the `geos`
/// feature returns -1 with no side effects.
///
/// # Safety
/// All pointer arguments must be valid for reads of their declared lengths
/// for the duration of the call: coord arrays hold coords_len doubles with
/// an even length of at least 6, hole arrays hold holes_len entries, and
/// out-pointers must be non-null with room for the documented result
/// struct. Null or short inputs return -1 instead of crashing.
#[no_mangle]
pub unsafe extern "C" fn mic_solve(
    coords: *const c_double,
    coords_len: size_t,
    holes: *const HoleInput,
    holes_len: size_t,
    options: *const MicOptionsC,
    result: *mut MicResultC,
) -> c_int {
    guard(|| {
        if result.is_null() {
            return -1;
        }
        let poly = match decode_polygon(coords, coords_len, holes, holes_len) {
            Some(p) => p,
            None => return -1,
        };
        let raw = if options.is_null() {
            MicOptionsC::default()
        } else {
            *options
        };
        #[cfg(not(feature = "geos"))]
        if raw.engine == 1 {
            return -1;
        }
        match maximum_inscribed_circle(&poly, &opts_from(raw)) {
            Ok(r) => fill(&r, result),
            Err(_) => -1,
        }
    })
}

/// Solve MIC over a multipolygon, keep the best component.
///
/// Returns 0 on success, -1 on error.
///
/// # Safety
/// Same contract as [`mic_solve`]: `polys` must point to `polys_len`
/// readable component structs, each with its own valid arrays.
#[no_mangle]
pub unsafe extern "C" fn mic_solve_multipolygon(
    polys: *const MicPolyInput,
    polys_len: size_t,
    options: *const MicOptionsC,
    result: *mut MicResultC,
) -> c_int {
    guard(|| {
        if polys.is_null() || result.is_null() {
            return -1;
        }
        let s = std::slice::from_raw_parts(polys, polys_len);
        let mut vec = Vec::with_capacity(s.len());
        for p in s {
            let poly = match decode_polygon(p.coords, p.coords_len, p.holes, p.holes_len) {
                Some(v) => v,
                None => return -1,
            };
            vec.push(poly);
        }
        let multi = geo_types::MultiPolygon(vec);
        let raw = if options.is_null() {
            MicOptionsC::default()
        } else {
            *options
        };
        #[cfg(not(feature = "geos"))]
        if raw.engine == 1 {
            return -1;
        }
        match maximum_inscribed_circle_multipolygon(&multi, &opts_from(raw)) {
            Ok(r) => fill(&r, result),
            Err(_) => -1,
        }
    })
}

#[no_mangle]
/// Return the default C ABI options: exact-then-GEOS engine, filtered mode.
pub extern "C" fn mic_options_default() -> MicOptionsC {
    MicOptionsC::default()
}
