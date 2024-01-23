//! corresponds to defs.h

/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h 20398 2016-11-22 17:17:12Z karypis $
 *
 */

use crate::{idx_t, real_t};

// pub const METISTITLE: &'static str ="METIS 5.2.1 Copyright 1998-22, Regents of the University of Minnesota\n";
// pub const MAXLINE: usize = 			1280000;

/// Initial number of maximum number of adjacent domains. This number will be adjusted as required.
// pub const INIT_MAXNAD: idx_t =              200;

// #define LARGENIPARTS		7	/* Number of random initial partitions */
// #define SMALLNIPARTS		5	/* Number of random initial partitions */
pub const COARSEN_FRACTION: real_t = 0.85; /* Node reduction between successive coarsening levels */

pub const HTLENGTH: idx_t = (1 << 13) - 1;
// pub const COMPRESSION_FRACTION: real_t = 		0.85;

// #define MMDSWITCH		        120

/* Default ufactors for the various operational modes */
// #define PMETIS_DEFAULT_UFACTOR          1
// #define MCPMETIS_DEFAULT_UFACTOR        10
// #define KMETIS_DEFAULT_UFACTOR          30
// #define OMETIS_DEFAULT_UFACTOR          200
