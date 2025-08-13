/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * timing.c
 *
 * This file contains routines that deal with timing Metis
 *
 * Started 7/24/97
 * George
 *
 * $Id: timing.c 13936 2013-03-30 03:59:09Z karypis $
 *
 */

use crate::*;

pub fn gk_clearcputimer(timer: &mut f64) {
    *timer = 0.0;
}

/// even though this is a nop, I want to keep all it's uses so I can deprecate it and replace them
pub fn gk_getcputimer(timer: f64) -> f64 {
    timer
}

// TODO: replace all instances with std::time::Duration
fn get_cpu_seconds() -> f64 {
    #[cfg(not(target_family = "unix"))]
    compile_error!("TODO: support windows time");

    #[cfg(target_family = "unix")]
    {
        let mut r: libc::rusage = unsafe { std::mem::zeroed() };
        unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut r) };
        (r.ru_utime.tv_sec + r.ru_stime.tv_sec) as f64  + 1.0e-6*(r.ru_utime.tv_usec + r.ru_stime.tv_usec) as f64
    }
}

pub fn gk_startcputimer(timer: &mut f64) {
    *timer -= get_cpu_seconds()
}

pub fn gk_stopcputimer(timer: &mut f64) {
    *timer += get_cpu_seconds()
}

/*************************************************************************
* This function clears the timers
**************************************************************************/
#[metis_func]
pub extern "C" fn InitTimers(ctrl: &mut ctrl_t) {
    gk_clearcputimer(&mut ctrl.TotalTmr);
    gk_clearcputimer(&mut ctrl.InitPartTmr);
    gk_clearcputimer(&mut ctrl.MatchTmr);
    gk_clearcputimer(&mut ctrl.ContractTmr);
    gk_clearcputimer(&mut ctrl.CoarsenTmr);
    gk_clearcputimer(&mut ctrl.UncoarsenTmr);
    gk_clearcputimer(&mut ctrl.RefTmr);
    gk_clearcputimer(&mut ctrl.ProjectTmr);
    gk_clearcputimer(&mut ctrl.SplitTmr);
    gk_clearcputimer(&mut ctrl.Aux1Tmr);
    gk_clearcputimer(&mut ctrl.Aux2Tmr);
    gk_clearcputimer(&mut ctrl.Aux3Tmr);
}

/*************************************************************************
* This function prints the various timers
**************************************************************************/
#[metis_func]
pub extern "C" fn PrintTimers(ctrl: &mut ctrl_t) {
    print!("\nTiming Information -------------------------------------------------");
    print!("\n Multilevel: \t\t {:7.3}", gk_getcputimer(ctrl.TotalTmr));
    print!(
        "\n     Coarsening: \t\t {:7.3}",
        gk_getcputimer(ctrl.CoarsenTmr)
    );
    print!(
        "\n            Matching: \t\t\t {:7.3}",
        gk_getcputimer(ctrl.MatchTmr)
    );
    print!(
        "\n            Contract: \t\t\t {:7.3}",
        gk_getcputimer(ctrl.ContractTmr)
    );
    print!(
        "\n     Initial Partition: \t {:7.3}",
        gk_getcputimer(ctrl.InitPartTmr)
    );
    print!(
        "\n     Uncoarsening: \t\t {:7.3}",
        gk_getcputimer(ctrl.UncoarsenTmr)
    );
    print!(
        "\n          Refinement: \t\t\t {:7.3}",
        gk_getcputimer(ctrl.RefTmr)
    );
    print!(
        "\n          Projection: \t\t\t {:7.3}",
        gk_getcputimer(ctrl.ProjectTmr)
    );
    print!(
        "\n     Splitting: \t\t {:7.3}",
        gk_getcputimer(ctrl.SplitTmr)
    );
    /*
      print!("\n       Aux1Tmr: \t\t {:7.3}", gk_getcputimer(ctrl.Aux1Tmr));
      print!("\n       Aux2Tmr: \t\t {:7.3}", gk_getcputimer(ctrl.Aux2Tmr));
      print!("\n       Aux3Tmr: \t\t {:7.3}", gk_getcputimer(ctrl.Aux3Tmr));
    */
    print!("\n********************************************************************\n");
}
