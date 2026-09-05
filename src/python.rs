//! Python bindings for the MIC solvers (feature `python`).
//!
//! Thin facade over the Rust API: exterior ring plus optional holes in,
//! result object out. Panics become `ValueError`, never a host crash.
//! The GEOS fallback stays compiled out unless the `geos` feature is
//! enabled; the engine param is validated with a clear error naming the
//! disabled path. Ported from IGE `ige-py` with paths pointed at this crate.

use super::mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon, MicEngine, MicOptions,
    MicUsedEngine, RobustMode,
};
use geo_types::{Coord, LineString, Polygon};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// One polygon inside a multipolygon input: exterior ring plus holes.
type MicPolygonInput = (Vec<(f64, f64)>, Option<Vec<Vec<(f64, f64)>>>);

/// Build a polygon from exterior ring plus optional holes.
pub fn build_polygon(
    exterior: Vec<(f64, f64)>,
    holes: Option<Vec<Vec<(f64, f64)>>>,
) -> PyResult<Polygon<f64>> {
    if exterior.len() < 3 {
        return Err(PyValueError::new_err(
            "polygon exterior must contain at least 3 points",
        ));
    }
    let coords: Vec<Coord<f64>> = exterior.into_iter().map(|(x, y)| Coord { x, y }).collect();
    let exterior_ls = LineString::from(coords);
    let interiors: Vec<LineString<f64>> = holes
        .unwrap_or_default()
        .into_iter()
        .map(|ring| {
            let ring_coords: Vec<Coord<f64>> =
                ring.into_iter().map(|(x, y)| Coord { x, y }).collect();
            LineString::from(ring_coords)
        })
        .collect();
    Ok(Polygon::new(exterior_ls, interiors))
}

/// Map a core error into a Python ValueError preserving the message.
pub fn core_err<E: std::fmt::Display>(e: E) -> PyErr {
    PyValueError::new_err(format!("solve failed: {e}"))
}

/// Run a closure behind the panic boundary. A Rust panic becomes a
/// Python ValueError instead of aborting the host process.
/// NOTE: effective only when the cdylib is built with panic=unwind.
pub fn guard<T>(f: impl FnOnce() -> PyResult<T> + std::panic::UnwindSafe) -> PyResult<T> {
    match std::panic::catch_unwind(f) {
        Ok(r) => r,
        Err(_) => Err(PyValueError::new_err(
            "solver panicked; input rejected without crashing host",
        )),
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyMicResult {
    #[pyo3(get)]
    pub center_x: f64,
    #[pyo3(get)]
    pub center_y: f64,
    #[pyo3(get)]
    pub radius: f64,
    #[pyo3(get)]
    pub radius_sq: f64,
    #[pyo3(get)]
    pub used_engine: String,
    #[pyo3(get)]
    pub candidate_count: usize,
    #[pyo3(get)]
    pub component_index: Option<usize>,
}

fn parse_engine(value: Option<&str>) -> PyResult<MicEngine> {
    match value.unwrap_or("exact_then_geos") {
        "exact_only" => Ok(MicEngine::ExactOnly),
        "fallback_only" => Ok(MicEngine::FallbackOnly),
        "exact_then_geos" => Ok(MicEngine::ExactThenGeos),
        other => Err(PyValueError::new_err(format!(
            "invalid engine '{other}'; expected exact_only|fallback_only|exact_then_geos"
        ))),
    }
}

fn parse_robust(value: Option<&str>) -> PyResult<RobustMode> {
    match value.unwrap_or("filtered") {
        "fast_f64" => Ok(RobustMode::FastF64),
        "filtered" => Ok(RobustMode::Filtered),
        other => Err(PyValueError::new_err(format!(
            "invalid robust_mode '{other}'; expected fast_f64|filtered"
        ))),
    }
}

fn engine_name(v: MicUsedEngine) -> &'static str {
    match v {
        MicUsedEngine::Exact => "exact",
        MicUsedEngine::GeosFallback => "geos_fallback",
        MicUsedEngine::Grid => "grid",
    }
}

fn check_geos(engine: MicEngine) -> PyResult<()> {
    #[cfg(not(feature = "geos"))]
    {
        if matches!(engine, MicEngine::FallbackOnly) {
            return Err(PyValueError::new_err(
                "engine 'fallback_only' needs the GEOS fallback, which is disabled in this wheel (built without the `geos` feature)",
            ));
        }
    }
    Ok(())
}

#[pyfunction(signature = (exterior, holes=None, engine=None, robust_mode=None))]
pub fn solve_mic_py(
    exterior: Vec<(f64, f64)>,
    holes: Option<Vec<Vec<(f64, f64)>>>,
    engine: Option<&str>,
    robust_mode: Option<&str>,
) -> PyResult<PyMicResult> {
    guard(|| {
        let engine_v = parse_engine(engine)?;
        check_geos(engine_v)?;
        let polygon = build_polygon(exterior, holes)?;
        let opts = MicOptions {
            engine: engine_v,
            use_bvh: false,
            use_parallel: false,
            robust_mode: parse_robust(robust_mode)?,
        };
        let r = maximum_inscribed_circle(&polygon, &opts).map_err(core_err)?;
        Ok(PyMicResult {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            used_engine: engine_name(r.used_engine).to_string(),
            candidate_count: r.candidate_count,
            component_index: r.component_index,
        })
    })
}

#[pyfunction(signature = (polygons, engine=None, robust_mode=None))]
pub fn solve_mic_multipolygon_py(
    polygons: Vec<MicPolygonInput>,
    engine: Option<&str>,
    robust_mode: Option<&str>,
) -> PyResult<PyMicResult> {
    guard(|| {
        let engine_v = parse_engine(engine)?;
        check_geos(engine_v)?;
        let robust = parse_robust(robust_mode)?;
        let mut polys = Vec::with_capacity(polygons.len());
        for (ext, holes) in polygons {
            polys.push(build_polygon(ext, holes)?);
        }
        let multi = geo_types::MultiPolygon(polys);
        let opts = MicOptions {
            engine: engine_v,
            use_bvh: false,
            use_parallel: false,
            robust_mode: robust,
        };
        let r = maximum_inscribed_circle_multipolygon(&multi, &opts).map_err(core_err)?;
        Ok(PyMicResult {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            used_engine: engine_name(r.used_engine).to_string(),
            candidate_count: r.candidate_count,
            component_index: r.component_index,
        })
    })
}

pub fn register(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let _ = py;
    m.add_class::<PyMicResult>()?;
    m.add_function(wrap_pyfunction!(solve_mic_py, m)?)?;
    m.add_function(wrap_pyfunction!(solve_mic_multipolygon_py, m)?)?;
    Ok(())
}
