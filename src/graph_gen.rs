use crate::*;

pub struct GraphBuilder {
    op: Optype,
    nparts: u32,
    ncon: u32,
    xadj: Vec<idx_t>,
    adjncy: Vec<idx_t>,
    adjwgt: Option<Vec<idx_t>>,
    vsize: Option<Vec<idx_t>>,
    vwgt: Option<Vec<idx_t>>,
    ubvec: Option<Vec<real_t>>,
    tpwgts: Option<Vec<real_t>>,
}

fn vec_ptr<T>(v: &mut Option<Vec<T>>) -> *mut T {
    match v.as_mut() {
        Some(v) => v.as_mut_ptr(),
        None => std::ptr::null_mut(),
    }
}

impl GraphBuilder {
    fn new(op: Optype, nparts: u32, ncon: u32) -> Self {
        Self {
            op,
            nparts,
            ncon,
            xadj: Vec::new(),
            adjncy: Vec::new(),
            adjwgt: None,
            vsize: None,
            vwgt: None,
            ubvec: None,
            tpwgts: None,
        }
    }

    /// Allocate every vertex to have corresponding vertex degree
    pub fn with_vtx_degrees(&mut self, mut deg: Vec<idx_t>) {
        debug_assert!(deg.iter().all(|&x| x >= 0), "cannot have negative degree");
        deg.push(0);
        util::make_csr(deg.len() - 1, &mut deg);
        self.xadj = deg;
    }

    /// using vertex degrees specified by `with_vtx_degrees`, make random edges
    pub fn random_edges(&mut self) {
        let mut first_open = 0;
    }
}
