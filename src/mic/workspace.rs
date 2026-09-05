use super::index::{NearestBoundaryIndex, PipIndex};
use super::input::{HostPolygon, SegmentIndex};
use super::MicError;

#[derive(Debug, Clone)]
pub(crate) struct MicCandidate {
    pub x: f64,
    pub y: f64,
    pub radius_sq: f64,
}

/// Reusable solver workspace for MIC computation.
#[derive(Debug, Clone)]
pub struct MicWorkspace {
    pub host: HostPolygon,
    pub pip_index: PipIndex,
    pub nb_index: NearestBoundaryIndex,
    pub(crate) candidate_buf: Vec<MicCandidate>,
}

impl MicWorkspace {
    pub fn new(host: HostPolygon, use_bvh: bool) -> Result<Self, MicError> {
        let seg_index = SegmentIndex::from_host(&host);
        if seg_index.is_empty() {
            return Err(MicError::InvalidInput(
                "polygon has no non-degenerate boundary segments".to_string(),
            ));
        }
        let pip_index = PipIndex::new(&host);
        let nb_index = NearestBoundaryIndex::new(seg_index, use_bvh);
        Ok(Self {
            host,
            pip_index,
            nb_index,
            candidate_buf: Vec::new(),
        })
    }

    pub(crate) fn clear_candidates(&mut self) {
        self.candidate_buf.clear();
    }
}
