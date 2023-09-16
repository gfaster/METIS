#![cfg(test)]
#![allow(unused_mut)]

use std::ptr;

use crate::bindings::{
    idx_t, METIS_PartGraphKway, METIS_PartGraphRecursive, METIS_NOPTIONS, METIS_OK,
};

use crate::util::{create_dummy_weights, verify_part};

#[test]
fn basic_part_graph_recursive() {
    let mut xadj = &mut [0, 1, 2];
    let mut adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let mut vwgt = ptr::null_mut();
    let mut vsize = ptr::null_mut();
    let mut adjwgt = ptr::null_mut();
    let mut nparts: idx_t = 2;
    let mut tpwgts = ptr::null_mut();
    let mut ubvec = ptr::null_mut();
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
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}

#[test]
fn basic_part_graph_kway() {
    let mut xadj = &mut [0, 1, 2];
    let mut adjncy = &mut [1, 0];
    let mut nvtxs = xadj.len() as idx_t - 1;
    let mut ncon: idx_t = 1;
    let mut vwgt = ptr::null_mut();
    let mut vsize = ptr::null_mut();
    let mut adjwgt = ptr::null_mut();
    let mut nparts: idx_t = 2;
    let mut tpwgts = ptr::null_mut();
    let mut ubvec = ptr::null_mut();
    let mut objval: idx_t = 0;
    let mut part = [0; 2];
    let mut options = [-1; METIS_NOPTIONS as usize];

    let res = unsafe {
        METIS_PartGraphKway(
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
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);
    assert_ne!(part[0], part[1]);
}

fn part_graph_and_verify(
    ncon: idx_t,
    options: &mut [i32; METIS_NOPTIONS as usize],
    nparts: idx_t,
    use_vwgt: bool,
    use_adjwgt: bool,
) {
    let (mut xadj, mut adjncy) = crate::util::read_graph(&mut std::io::Cursor::new(
        include_bytes!("../graphs/4elt_rs.graph"),
    ))
    .unwrap();

    let (mut v_vwgt, mut v_adjwgt) = create_dummy_weights(ncon as usize, &xadj, &adjncy);
    let vwgt;
    match use_vwgt {
        true => vwgt = v_vwgt.as_mut_ptr(),
        false => vwgt = ptr::null_mut(),
    }
    let adjwgt;
    match use_adjwgt {
        true => adjwgt = v_adjwgt.as_mut_ptr(),
        false => adjwgt = ptr::null_mut(),
    }
    let mut ncon = ncon;

    let mut nparts = nparts;
    let mut part = vec![-1; xadj.len() - 1];

    let mut objval = 0;

    let res = unsafe {
        METIS_PartGraphKway(
            &mut (xadj.len() as idx_t - 1) as *mut idx_t,
            &mut ncon as *mut _,
            xadj.as_mut_ptr(),
            adjncy.as_mut_ptr(),
            vwgt,
            ptr::null_mut(),
            adjwgt,
            &mut nparts as *mut _,
            ptr::null_mut(),
            ptr::null_mut(),
            options.as_mut_ptr(),
            &mut objval as *mut _,
            part.as_mut_ptr(),
        )
    };

    assert_eq!(res, METIS_OK);

    let vwgt = match use_vwgt {
        true => Some(&v_vwgt[..]),
        false => None,
    };
    let adjwgt = match use_adjwgt {
        true => Some(&v_adjwgt[..]),
        false => None,
    };

    verify_part(&xadj, &adjncy, vwgt, adjwgt, objval, &part, nparts)
}

macro_rules! part_test {
    (
    name: $name:ident,
    options: $options:expr,
    nparts: $nparts:literal,
    ncon: $ncon:literal,
    vwgt: $use_vwgt:literal,
    adjwgt: $use_adjwgt:literal,
    ) => {
        #[test]
        fn $name() {
            let mut options = $options;
            part_graph_and_verify($ncon, &mut options, $nparts, $use_vwgt, $use_adjwgt);
        }
    };
}

part_test! {
    name: large_kway_1con_vwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 1,
    vwgt: true,
    adjwgt: false,
}

part_test! {
    name: large_kway_2con_vwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 2,
    vwgt: true,
    adjwgt: false,
}

part_test! {
    name: large_kway_1con_vwgt_adjwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 1,
    vwgt: true,
    adjwgt: true,
}

part_test! {
    name: large_kway_2con_vwgt_adjwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 2,
    vwgt: true,
    adjwgt: true,
}

part_test! {
    name: large_kway_halve,
    options: make_options!(Cut Grow),
    nparts: 2,
    ncon: 1,
    vwgt: false,
    adjwgt: false,
}

part_test! {
    name: large_kway_1con_adjwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 1,
    vwgt: false,
    adjwgt: true,
}

part_test! {
    name: large_kway_2con_adjwgt,
    options: make_options!(Cut Grow),
    nparts: 20,
    ncon: 2,
    vwgt: false,
    adjwgt: true,
}

part_test! {
    name: large_kway_vol_basic,
    options: make_options!(Vol Grow),
    nparts: 20,
    ncon: 1,
    vwgt: false,
    adjwgt: false,
}

part_test! {
    name: large_kway_vol_vwgt,
    options: make_options!(Vol Grow),
    nparts: 20,
    ncon: 1,
    vwgt: true,
    adjwgt: false,
}

part_test! {
    name: large_kway_vol_2con_vwgt,
    options: make_options!(Vol Grow),
    nparts: 20,
    ncon: 2,
    vwgt: true,
    adjwgt: false,
}
