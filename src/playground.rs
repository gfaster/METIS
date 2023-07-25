#[allow(unused_imports)]
use metis::*;

fn main() {
    unimplemented!();
}



#[test]
fn playground_can_utilize_lib() {
    let xadj = &mut [0, 1, 2];
    let adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let vwgt = std::ptr::null_mut();
    let vsize = std::ptr::null_mut();
    let adjwgt = std::ptr::null_mut();
    let mut nparts: idx_t = 2;
    let tpwgts = std::ptr::null_mut();
    let ubvec = std::ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        METIS_PartGraphRecursive(
            &mut nvtxs as *mut _,
            &mut ncon as *mut _,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            vsize,
            adjwgt,
            &mut nparts as *mut _,
            tpwgts,
            ubvec,
            options.as_mut_ptr(),
            &mut objval as *mut _,
            part.as_mut_ptr()
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}
