//! MIC tuning constants carried over from IGE `ige-core::tuning`.
//!
//! Only the `MIC_*` family. No grid, expand, certify, or nesting constants.

pub const MIC_BINARY_ITERS: usize = 40;
pub const MIC_EXPANSION_ITERS: usize = 50;
pub const MIC_EPS: f64 = 1e-12;
pub const MIC_CANDIDATE_QUANTIZE: f64 = 1e9;
pub const MIC_BASE_TRIPLE_CAP: usize = 64;
pub const MIC_BASE_SS_SEG_CAP: usize = 32;
pub const MIC_BASE_SS_VERT_CAP: usize = 12;
pub const MIC_BASE_SEGS_PER_RING: usize = 3;
pub const MIC_EXT_TRIPLE_CAP: usize = 96;
pub const MIC_EXT_SS_SEG_CAP: usize = 64;
pub const MIC_EXT_SS_VERT_CAP: usize = 32;
pub const MIC_EXT_SEGS_PER_RING: usize = 5;
