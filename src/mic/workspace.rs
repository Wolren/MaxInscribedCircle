//! Reusable solver workspace that owns the host polygon and the indexes
//! built for it.
use super::index::{NearestBoundaryIndex, PipIndex};
use super::input::{HostPolygon, SegmentIndex};
use super::MicError;

/// Reusable solver workspace for MIC computation.
#[derive(Debug, Clone)]
pub struct MicWorkspace {
    /// Host polygon with its flattened ring metadata.
    pub host: HostPolygon,
    /// Point-in-polygon index over the host rings.
    pub pip_index: PipIndex,
    /// Nearest-boundary distance index over the host segments.
    pub nb_index: NearestBoundaryIndex,
}

impl MicWorkspace {
    /// Build both indexes for `host`.
    pub fn new(host: HostPolygon) -> Result<Self, MicError> {
        let seg_index = SegmentIndex::from_host(&host);
        Self::with_segments(host, seg_index)
    }

    /// Build both indexes from a segment table the caller already has.
    pub fn with_segments(host: HostPolygon, seg_index: SegmentIndex) -> Result<Self, MicError> {
        if seg_index.is_empty() {
            return Err(MicError::InvalidInput(
                "polygon has no non-degenerate boundary segments".to_string(),
            ));
        }
        let pip_index = PipIndex::new(&host);
        let nb_index = NearestBoundaryIndex::new(seg_index);
        Ok(Self {
            host,
            pip_index,
            nb_index,
        })
    }
}
