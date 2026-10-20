//! Python bindings for the MIC solvers (feature `python`).
//!
//! Thin facade over the Rust API: exterior ring plus optional holes in,
//! result object out. Panics become `ValueError`, never a host crash.
//! Ported from IGE `ige-py` with paths pointed at this crate.

use super::mic::{
    maximum_inscribed_circle, maximum_inscribed_circle_multipolygon, MicOptions, MicUsedEngine,
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

#[pyclass(skip_from_py_object)]
#[derive(Clone)]
/// MIC solve result exposed to Python callers.
pub struct PyMicResult {
    #[pyo3(get)]
    /// Circle center x coordinate.
    pub center_x: f64,
    #[pyo3(get)]
    /// Circle center y coordinate.
    pub center_y: f64,
    #[pyo3(get)]
    /// Circle radius in input coordinate units.
    pub radius: f64,
    #[pyo3(get)]
    /// Squared radius, radius * radius.
    pub radius_sq: f64,
    #[pyo3(get)]
    /// Solver that produced the result: exact or grid.
    pub used_engine: String,
    #[pyo3(get)]
    /// Number of candidate centers evaluated during the solve.
    pub candidate_count: usize,
    #[pyo3(get)]
    /// Index of the multipolygon component that holds the circle, if any.
    pub component_index: Option<usize>,
}

fn used_engine_name(v: MicUsedEngine) -> &'static str {
    match v {
        MicUsedEngine::Grid => "grid",
        MicUsedEngine::Geos => "geos",
    }
}

#[pyfunction(signature = (exterior, holes=None))]
/// Solve MIC for one polygon given as an exterior ring plus optional holes.
pub fn solve_mic_py(
    exterior: Vec<(f64, f64)>,
    holes: Option<Vec<Vec<(f64, f64)>>>,
) -> PyResult<PyMicResult> {
    guard(|| {
        let polygon = build_polygon(exterior, holes)?;
        let opts = MicOptions::default();
        let r = maximum_inscribed_circle(&polygon, &opts).map_err(core_err)?;
        Ok(PyMicResult {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            used_engine: used_engine_name(r.used_engine).to_string(),
            candidate_count: r.candidate_count,
            component_index: r.component_index,
        })
    })
}

#[pyfunction(signature = (polygons))]
/// Solve MIC across several polygons and return the best component circle.
pub fn solve_mic_multipolygon_py(polygons: Vec<MicPolygonInput>) -> PyResult<PyMicResult> {
    guard(|| {
        let mut polys = Vec::with_capacity(polygons.len());
        for (ext, holes) in polygons {
            polys.push(build_polygon(ext, holes)?);
        }
        let multi = geo_types::MultiPolygon(polys);
        let opts = MicOptions::default();
        let r = maximum_inscribed_circle_multipolygon(&multi, &opts).map_err(core_err)?;
        Ok(PyMicResult {
            center_x: r.center.x(),
            center_y: r.center.y(),
            radius: r.radius,
            radius_sq: r.radius_sq,
            used_engine: used_engine_name(r.used_engine).to_string(),
            candidate_count: r.candidate_count,
            component_index: r.component_index,
        })
    })
}

/// Register the result class and the solver functions with the module.
pub fn register(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    let _ = py;
    m.add_class::<PyMicResult>()?;
    m.add_function(wrap_pyfunction!(solve_mic_py, m)?)?;
    m.add_function(wrap_pyfunction!(solve_mic_multipolygon_py, m)?)?;
    Ok(())
}
