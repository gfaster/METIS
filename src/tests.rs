#![cfg(test)]
#![allow(unused_mut)]

use std::ptr;

use crate::bindings::{idx_t, real_t, METIS_NOPTIONS, METIS_OK};
use crate::kmetis::METIS_PartGraphKway;
use crate::pmetis::METIS_PartGraphRecursive;

use crate::util::{create_dummy_weights, verify_part};

/// function signature of METIS_PartGraphKway, METIS_PartGraphRecursive
type PartSig = unsafe extern "C" fn(
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
    *mut real_t,
    *const real_t,
    *mut idx_t,
    *mut idx_t,
    *mut idx_t,
) -> idx_t;

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
    partfn: PartSig,
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
        partfn(
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
    partfn: $partfn:ident,
    extra: $extra:ident,
    ) => {
        #[cfg_attr(not(feature = "extra_tests"), ignore = "extra test")]
        #[test]
        fn $name() {
            let mut options = $options;
            part_graph_and_verify(
                $ncon,
                &mut options,
                $nparts,
                $use_vwgt,
                $use_adjwgt,
                $partfn,
            );
        }
    };
    (
    name: $name:ident,
    options: $options:expr,
    nparts: $nparts:literal,
    ncon: $ncon:literal,
    vwgt: $use_vwgt:literal,
    adjwgt: $use_adjwgt:literal,
    partfn: $partfn:ident,
    ) => {
        #[test]
        fn $name() {
            let mut options = $options;
            part_graph_and_verify(
                $ncon,
                &mut options,
                $nparts,
                $use_vwgt,
                $use_adjwgt,
                $partfn,
            );
        }
    };
}

macro_rules! part_test_set {
    (
    set_name: $name:ident,
    options: $options:expr,
    partfn: $partfn:ident,
    ) => {
        mod $name {
            use super::*;

            part_test! {
                name: large_basic,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: false,
                adjwgt: false,
                partfn: $partfn,
            }

            part_test! {
                name: large_1con_vwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: true,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_2con_vwgt,
                options: $options,
                nparts: 20,
                ncon: 2,
                vwgt: true,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_1con_vwgt_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: true,
                adjwgt: true,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_2con_vwgt_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 2,
                vwgt: true,
                adjwgt: true,
                partfn: $partfn,
            }

            part_test! {
                name: large_halve,
                options: $options,
                nparts: 2,
                ncon: 1,
                vwgt: false,
                adjwgt: false,
                partfn: $partfn,
                extra: true,
            }

            part_test! {
                name: large_1con_adjwgt,
                options: $options,
                nparts: 20,
                ncon: 1,
                vwgt: false,
                adjwgt: true,
                partfn: $partfn,
                extra: true,
            }
        }
    };
}

// this could be shrunk but it works
macro_rules! part_test_hyper_set {
    ($partfn:ident as $name:ident => [$({$($first_opt:ident),*}),*]) => {
        mod $name {
            #![allow(non_snake_case, unused_imports)]
            use super::*;
            part_test_hyper_set!(@call_for $partfn : ($(($($first_opt),*)),*) @ ());
        }
    };
    (@call $partfn:ident, $name:ident: ($($option:ident),*)) => {
        // $(part_test_hyper_set!(@dbg $option));*
        part_test_set!(
            set_name: $name,
            options: make_options!($($option)*),
            partfn: $partfn,
        );
    };
    (@call_for $partfn:ident: ($first:tt $(, $rest:tt)*) @ $keep:tt) => {
                // $(part_test_hyper_set!(@dbg $rest))*;
                // part_test_hyper_set!(@dbg $keep);
        part_test_hyper_set!(@call_for_inner $partfn: $first @ ($($rest),*) @ $keep);
    };
    (@call_for_inner $partfn:ident: ($($iter:ident),*) @ () @ $keep:tt) => {
        $(
        // mod $iter {
        //     use super::*;
        // part_test_hyper_set!(@dbg $iter);
            part_test_hyper_set!(@rejoin_call $partfn, $iter: $keep);
        // }
        )*
    };
    (@call_for_inner $partfn:ident: ($($iter:ident),*) @ $rest:tt @ $keep:tt) => {
        $(
            mod $iter {
                use super::*;
        // part_test_hyper_set!(@dbg $keep);
                part_test_hyper_set!(@rejoin $partfn: $rest @ ($iter, $keep));
            }
        )*
    };
    (@rejoin $partfn:ident: $rest:tt @ ($parent:ident, ($($keep:ident),*))) => {
        // part_test_hyper_set!(@dbg $parent);
        // part_test_hyper_set!(@dbg $rest);
        part_test_hyper_set!(@call_for $partfn: $rest @ ($parent $(, $keep)*));
    };
    (@rejoin_call $partfn:ident, $name:ident: ($($keep:ident),*)) => {
        part_test_hyper_set!(@call $partfn, $name: ($name $(, $keep)*));
    };
    (@dbg $var:expr) => {
        part_test_hyper_set!(@dbg_inner $var);
    };
    (@dbg_inner a) => {
    };
}

// TODO: test graphs with odd degree (very high, power law)
// I suspect problems if there is vertex degree > nparts

// Volume + Contig fails even when my code isn't used at all. I should verify this later, but it
// appears to be an issue in the original as well.

// Kmetis only allows Grow and Rb initial partitioning
part_test_hyper_set!(METIS_PartGraphKway as kmetisvol => [{Vol}, {Grow, Rb}, {Rm, Shem}]);
part_test_hyper_set!(METIS_PartGraphKway as kmetiscut => [{Cut}, {Grow, Rb}, {Rm, Shem}, {Contig, None}]);

// Communication volume is illegal on PMETIS routines
// Edge and Node initial partitioning also illegal in PMETIS
part_test_hyper_set!(METIS_PartGraphRecursive as pmetis => [{Cut}, {Grow, Random}, {Rm, Shem}]);
