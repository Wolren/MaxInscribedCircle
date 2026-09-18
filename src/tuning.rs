//! MIC tuning constants carried over from IGE `ige-core::tuning`.
//!
//! Only the `MIC_*` family. No grid, expand, certify, or nesting constants.

/// Binary search iterations for the convex quadrilateral Chebyshev radius.
pub const MIC_BINARY_ITERS: usize = 40;
/// Binary search iterations for the convex n-gon feasible radius.
pub const MIC_EXPANSION_ITERS: usize = 50;
/// Determinant magnitude below which a linear solve counts as singular.
pub const MIC_EPS: f64 = 1e-12;
/// Candidate deduplication scale: coordinates are multiplied by this and rounded.
pub const MIC_CANDIDATE_QUANTIZE: f64 = 1e9;
/// Maximum segments sampled for segment-triple candidates on simple polygons.
pub const MIC_BASE_TRIPLE_CAP: usize = 64;
/// Maximum segments sampled for seg-seg-vertex candidates on simple polygons.
pub const MIC_BASE_SS_SEG_CAP: usize = 32;
/// Maximum vertices sampled for seg-seg-vertex candidates on simple polygons.
pub const MIC_BASE_SS_VERT_CAP: usize = 12;
/// Per-ring segment budget of the base sampling tier, ignored by the sampler.
pub const MIC_BASE_SEGS_PER_RING: usize = 3;
/// Maximum segments sampled for segment-triple candidates on complex polygons.
pub const MIC_EXT_TRIPLE_CAP: usize = 96;
/// Maximum segments sampled for seg-seg-vertex candidates on complex polygons.
pub const MIC_EXT_SS_SEG_CAP: usize = 64;
/// Maximum vertices sampled for seg-seg-vertex candidates on complex polygons.
pub const MIC_EXT_SS_VERT_CAP: usize = 32;
/// Per-ring segment budget of the extended tier, ignored by the sampler.
pub const MIC_EXT_SEGS_PER_RING: usize = 5;
/// Half of the first-level cell count of the grid solver: the bounding box is
/// covered by at most this many cells along its longer axis.
pub const MIC_GRID_SIDE: usize = 25;
/// Relative gap above which a fast-path radius counts as unverified. A radius
/// the segment table does not confirm to this accuracy means the polygon is not
/// the class the fast path assumed, so the general path has to answer instead.
pub const MIC_FAST_PATH_REL_TOL: f64 = 1e-9;
