#![allow(bad_style, unused)]

#[metis_decl]
extern "C" {
    pub fn SetupCtrl(
        optype: moptype_et,
        options: *const idx_t, // manually verified const
        ncon: idx_t,
        nparts: idx_t,
        tpwgts: *const real_t, // manually verified const
        ubvec: *const real_t,  // manually verified const
    ) -> *mut ctrl_t;

    pub fn AllocateRefinementWorkSpace(
        ctrl: *mut ctrl_t,
        nbrpoolsize_max: idx_t,
        nbrpoolsize: idx_t,
    ) -> std::ffi::c_void;
    pub fn Allocate2WayNodePartitionMemory(ctrl: *mut ctrl_t, graph: *mut graph_t) -> std::ffi::c_void;

    pub fn AllocateWorkSpace(ctrl: *mut ctrl_t, graph: *mut graph_t) -> std::ffi::c_void;
    pub fn cnbrpoolGetNext(ctrl: *mut ctrl_t, nnbrs: idx_t) -> idx_t;
    pub fn cnbrpoolReset(ctrl: *mut ctrl_t) -> std::ffi::c_void;

    pub fn imalloc(nmemb: usize, msg: *const std::ffi::c_char) -> *mut std::ffi::c_void;
    pub fn ismalloc(
        nmemb: usize,
        val: idx_t,
        msg: *const std::ffi::c_char,
    ) -> *mut std::ffi::c_void;
    pub fn rsmalloc(
        nmemb: usize,
        val: real_t,
        msg: *const std::ffi::c_char,
    ) -> *mut std::ffi::c_void;
    pub fn rmalloc(nmemb: usize, msg: *const std::ffi::c_char) -> *mut std::ffi::c_void;
    pub fn irealloc(
        old: *mut std::ffi::c_void,
        nmemb: usize,
        msg: *const std::ffi::c_char,
    ) -> *mut std::ffi::c_void;
    pub fn FreeCtrl(r_ctrl: *mut *mut ctrl_t) -> std::ffi::c_void;
    pub fn FreeWorkSpace(ctrl: *mut ctrl_t) -> std::ffi::c_void;

    pub fn isrand(seed: idx_t) -> std::ffi::c_void;
    pub fn SetupKWayBalMultipliers(ctrl: *mut ctrl_t, graph: *mut graph_t) -> std::ffi::c_void;
    pub fn vnbrpoolGetNext(ctrl: *mut ctrl_t, nnbrs: idx_t) -> idx_t;
    pub fn vnbrpoolReset(ctrl: *mut ctrl_t) -> std::ffi::c_void;

    pub fn irandArrayPermute(
        n: idx_t,
        p: *mut idx_t,
        nshuffles: idx_t,
        flag: std::ffi::c_int,
    ) -> std::ffi::c_void;
    pub fn irandInRange(r: idx_t) -> idx_t;

    pub fn Setup2WayBalMultipliers(
        ctrl: *mut ctrl_t,
        graph: *mut graph_t,
        tpwgts: *mut real_t,
    ) -> std::ffi::c_void;
    pub fn Compute2WayNodePartitionParams(
        ctrl: *mut ctrl_t,
        graph: *mut graph_t,
    ) -> std::ffi::c_void;
    pub fn FM_2WayNodeRefine2Sided(
        ctrl: *mut ctrl_t,
        graph: *mut graph_t,
        niter: idx_t,
    ) -> std::ffi::c_void;
    pub fn FM_2WayNodeRefine1Sided(
        ctrl: *mut ctrl_t,
        graph: *mut graph_t,
        niter: idx_t,
    ) -> std::ffi::c_void;


    pub fn Refine2WayNode(ctrl: *mut ctrl_t, orginal_graph: *mut graph_t, graph: *mut graph_t);
}


// replacing c args with rust args (2 commands over visual region)
// 1:
// '<,'>s/\v(\w+_t) (\**)(\w+)/\3: \2\1/g
// 2:
// '<,'>s/\V*/*mut /g

pub const METIS_VER_MAJOR: u32 = 5;
pub const METIS_VER_MINOR: u32 = 1;
pub const METIS_VER_SUBMINOR: u32 = 0;
pub const METIS_NOPTIONS: u32 = 40;
pub type idx_t = i32;
pub type real_t = f32;
pub type ipq_t = std::ffi::c_void;

// gklib functions are not renamed
extern "C" {
    /// initialize malloc used by most METIS functions
    /// it sets a thread-local variable for the core, so it should be fine to use in tests
    /// any gk_malloc operations made after calling init will be freed, and calling gk_free on any
    /// gk_malloc calls made before are not valid to be freed.
    pub fn gk_malloc_init() -> std::ffi::c_int;
    pub fn gk_malloc_cleanup(showstats: std::ffi::c_int) -> std::ffi::c_void;
    pub fn gk_malloc(size: usize, msg: *const std::ffi::c_char) -> *mut std::ffi::c_void;

    /// my wrapper for gk_free that isn't variadic
    pub fn gk_free_one(ptr: *mut *mut std::ffi::c_void) -> std::ffi::c_void;
}

// these don't need metis_decl attrib since they are the public API
// extern "C" {
//     pub fn METIS_PartGraphRecursive(
//         nvtxs: *mut idx_t,
//         ncon: *mut idx_t,
//         xadj: *mut idx_t,
//         adjncy: *mut idx_t,
//         vwgt: *mut idx_t,
//         vsize: *mut idx_t,
//         adjwgt: *mut idx_t,
//         nparts: *mut idx_t,
//         tpwgts: *mut real_t,
//         ubvec: *const real_t, // only copied to internal ubfactors array
//         options: *mut idx_t,
//         edgecut: *mut idx_t,
//         part: *mut idx_t,
//     ) -> ::std::os::raw::c_int;
// }
// extern "C" {
//     pub fn METIS_PartGraphKway(
//         nvtxs: *mut idx_t,
//         ncon: *mut idx_t,
//         xadj: *mut idx_t,
//         adjncy: *mut idx_t,
//         vwgt: *mut idx_t,
//         vsize: *mut idx_t,
//         adjwgt: *mut idx_t,
//         nparts: *mut idx_t,
//         tpwgts: *mut real_t,
//         ubvec: *mut real_t,
//         options: *mut idx_t,
//         edgecut: *mut idx_t,
//         part: *mut idx_t,
//     ) -> ::std::os::raw::c_int;
// }
extern "C" {
    pub fn METIS_MeshToDual(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        ncommon: *mut idx_t,
        numflag: *mut idx_t,
        r_xadj: *mut *mut idx_t,
        r_adjncy: *mut *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_MeshToNodal(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        numflag: *mut idx_t,
        r_xadj: *mut *mut idx_t,
        r_adjncy: *mut *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_PartMeshNodal(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        options: *mut idx_t,
        objval: *mut idx_t,
        epart: *mut idx_t,
        npart: *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_PartMeshDual(
        ne: *mut idx_t,
        nn: *mut idx_t,
        eptr: *mut idx_t,
        eind: *mut idx_t,
        vwgt: *mut idx_t,
        vsize: *mut idx_t,
        ncommon: *mut idx_t,
        nparts: *mut idx_t,
        tpwgts: *mut real_t,
        options: *mut idx_t,
        objval: *mut idx_t,
        epart: *mut idx_t,
        npart: *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_Free(ptr: *mut ::std::os::raw::c_void) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_SetDefaultOptions(options: *mut idx_t) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_NodeNDP(
        nvtxs: idx_t,
        xadj: *mut idx_t,
        adjncy: *mut idx_t,
        vwgt: *mut idx_t,
        npes: idx_t,
        options: *mut idx_t,
        perm: *mut idx_t,
        iperm: *mut idx_t,
        sizes: *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_ComputeVertexSeparator(
        nvtxs: *mut idx_t,
        xadj: *mut idx_t,
        adjncy: *mut idx_t,
        vwgt: *mut idx_t,
        options: *mut idx_t,
        sepsize: *mut idx_t,
        part: *mut idx_t,
    ) -> ::std::os::raw::c_int;
}
extern "C" {
    pub fn METIS_NodeRefine(
        nvtxs: idx_t,
        xadj: *mut idx_t,
        vwgt: *mut idx_t,
        adjncy: *mut idx_t,
        where_: *mut idx_t,
        hmarker: *mut idx_t,
        ubfactor: real_t,
    ) -> ::std::os::raw::c_int;
}
/// SIGABRT
pub const SIGMEM: std::os::raw::c_int = 6;

pub const METIS_OK: rstatus_et = 1;
pub const METIS_ERROR_INPUT: rstatus_et = -2;
pub const METIS_ERROR_MEMORY: rstatus_et = -3;
pub const METIS_ERROR: rstatus_et = -4;
pub type rstatus_et = ::std::os::raw::c_int;
pub const METIS_OP_PMETIS: moptype_et = 0;
pub const METIS_OP_KMETIS: moptype_et = 1;
pub const METIS_OP_OMETIS: moptype_et = 2;
pub type moptype_et = ::std::os::raw::c_uint;
pub const METIS_OPTION_PTYPE: moptions_et = 0;
pub const METIS_OPTION_OBJTYPE: moptions_et = 1;
pub const METIS_OPTION_CTYPE: moptions_et = 2;
pub const METIS_OPTION_IPTYPE: moptions_et = 3;
pub const METIS_OPTION_RTYPE: moptions_et = 4;
pub const METIS_OPTION_DBGLVL: moptions_et = 5;
pub const METIS_OPTION_NIPARTS: moptions_et = 6;
pub const METIS_OPTION_NITER: moptions_et = 7;
pub const METIS_OPTION_NCUTS: moptions_et = 8;
pub const METIS_OPTION_SEED: moptions_et = 9;
pub const METIS_OPTION_ONDISK: moptions_et = 10;
pub const METIS_OPTION_MINCONN: moptions_et = 11;
pub const METIS_OPTION_CONTIG: moptions_et = 12;
pub const METIS_OPTION_COMPRESS: moptions_et = 13;
pub const METIS_OPTION_CCORDER: moptions_et = 14;
pub const METIS_OPTION_PFACTOR: moptions_et = 15;
pub const METIS_OPTION_NSEPS: moptions_et = 16;
pub const METIS_OPTION_UFACTOR: moptions_et = 17;
pub const METIS_OPTION_NUMBERING: moptions_et = 18;
pub const METIS_OPTION_DROPEDGES: moptions_et = 19;
pub const METIS_OPTION_NO2HOP: moptions_et = 20;
pub const METIS_OPTION_TWOHOP: moptions_et = 21;
pub const METIS_OPTION_FAST: moptions_et = 22;

pub const METIS_OPTION_HELP: moptions_et = 23;
pub const METIS_OPTION_TPWGTS: moptions_et = 24;
pub const METIS_OPTION_NCOMMON: moptions_et = 25;
pub const METIS_OPTION_NOOUTPUT: moptions_et = 26;
pub const METIS_OPTION_BALANCE: moptions_et = 27;
pub const METIS_OPTION_GTYPE: moptions_et = 28;
pub const METIS_OPTION_UBVEC: moptions_et = 29;

pub type moptions_et = ::std::os::raw::c_uint;
pub const METIS_PTYPE_RB: mptype_et = 0;
pub const METIS_PTYPE_KWAY: mptype_et = 1;
pub type mptype_et = ::std::os::raw::c_uint;
pub const METIS_GTYPE_DUAL: mgtype_et = 0;
pub const METIS_GTYPE_NODAL: mgtype_et = 1;
pub type mgtype_et = ::std::os::raw::c_uint;
pub const METIS_CTYPE_RM: mctype_et = 0;
pub const METIS_CTYPE_SHEM: mctype_et = 1;
pub type mctype_et = ::std::os::raw::c_uint;
pub const METIS_IPTYPE_GROW: miptype_et = 0;
pub const METIS_IPTYPE_RANDOM: miptype_et = 1;
pub const METIS_IPTYPE_EDGE: miptype_et = 2;
pub const METIS_IPTYPE_NODE: miptype_et = 3;
pub const METIS_IPTYPE_METISRB: miptype_et = 4;
pub type miptype_et = ::std::os::raw::c_uint;
pub const METIS_RTYPE_FM: mrtype_et = 0;
pub const METIS_RTYPE_GREEDY: mrtype_et = 1;
pub const METIS_RTYPE_SEP2SIDED: mrtype_et = 2;
pub const METIS_RTYPE_SEP1SIDED: mrtype_et = 3;
pub type mrtype_et = ::std::os::raw::c_uint;
pub const METIS_DBG_INFO: mdbglvl_et = 1;
pub const METIS_DBG_TIME: mdbglvl_et = 2;
pub const METIS_DBG_COARSEN: mdbglvl_et = 4;
pub const METIS_DBG_REFINE: mdbglvl_et = 8;
pub const METIS_DBG_IPART: mdbglvl_et = 16;
pub const METIS_DBG_MOVEINFO: mdbglvl_et = 32;
pub const METIS_DBG_SEPINFO: mdbglvl_et = 64;
pub const METIS_DBG_CONNINFO: mdbglvl_et = 128;
pub const METIS_DBG_CONTIGINFO: mdbglvl_et = 256;
pub const METIS_DBG_MEMORY: mdbglvl_et = 2048;
pub type mdbglvl_et = ::std::os::raw::c_uint;
pub const METIS_OBJTYPE_CUT: mobjtype_et = 0;
pub const METIS_OBJTYPE_VOL: mobjtype_et = 1;
pub const METIS_OBJTYPE_NODE: mobjtype_et = 2;
pub type mobjtype_et = ::std::os::raw::c_uint;

#[repr(u32)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DbgLvl {
    Info = METIS_DBG_INFO, 
    Time = METIS_DBG_TIME, 
    Coarsen = METIS_DBG_COARSEN, 
    Refine = METIS_DBG_REFINE, 
    Ipart = METIS_DBG_IPART, 
    MoveInfo = METIS_DBG_MOVEINFO, 
    SepInfo = METIS_DBG_SEPINFO, 
    ConnInfo = METIS_DBG_CONNINFO, 
    ContigInfo = METIS_DBG_CONTIGINFO, 
    Memory = METIS_DBG_MEMORY, 
}

// interals
pub(crate) const BNDTYPE_REFINE: std::ffi::c_int = 1;
pub(crate) const BNDTYPE_BALANCE: std::ffi::c_int = 2;

pub(crate) const OMODE_REFINE: std::ffi::c_int = 1;
pub(crate) const OMODE_BALANCE: std::ffi::c_int = 2;

pub(crate) const COMPRESSION_FRACTION: real_t = 0.85;

pub(crate) const LARGENIPARTS: idx_t = 7;
pub(crate) const SMALLNIPARTS: idx_t = 5;

/// unused except in binaries
#[repr(u32)]
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum Optype {
    Pmetis = METIS_OP_PMETIS,
    Kmetis = METIS_OP_KMETIS,
    Ometis = METIS_OP_OMETIS,
}

impl Optype {
    /// Returns `true` if the optype is [`Pmetis`].
    ///
    /// [`Pmetis`]: Optype::Pmetis
    #[must_use]
    pub fn is_pmetis(&self) -> bool {
        matches!(self, Self::Pmetis)
    }

    /// Returns `true` if the optype is [`Kmetis`].
    ///
    /// [`Kmetis`]: Optype::Kmetis
    #[must_use]
    pub fn is_kmetis(&self) -> bool {
        matches!(self, Self::Kmetis)
    }

    /// Returns `true` if the optype is [`Ometis`].
    ///
    /// [`Ometis`]: Optype::Ometis
    #[must_use]
    pub fn is_ometis(&self) -> bool {
        matches!(self, Self::Ometis)
    }
}

/// Fun fact: this mean "Objective Type", not "Object Type"
#[repr(u32)]
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum Objtype {
    Cut = METIS_OBJTYPE_CUT,
    Vol = METIS_OBJTYPE_VOL,
    // basically unused
    Node = METIS_OBJTYPE_NODE,
}

#[repr(u32)]
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum Iptype {
    Grow = METIS_IPTYPE_GROW,
    Random = METIS_IPTYPE_RANDOM,
    Edge = METIS_IPTYPE_EDGE,
    Node = METIS_IPTYPE_NODE,
    Rb = METIS_IPTYPE_METISRB,
}

/// only used in node recursive dissection
#[repr(u32)]
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum Rtype {
    Greedy = METIS_RTYPE_GREEDY,
    Fm = METIS_RTYPE_FM,
    Sep2Sided = METIS_RTYPE_SEP2SIDED,
    Sep1Sided = METIS_RTYPE_SEP1SIDED,
}

#[repr(u32)]
#[derive(PartialEq, Eq, Clone, Copy)]
pub enum Ctype {
    Rm = METIS_CTYPE_RM,
    Shem = METIS_CTYPE_SHEM,
}

/// we won't be using core
pub type gk_mcore_t = std::os::raw::c_void;

#[repr(C)]
pub struct ctrl_t {
    /// Type of operation
    pub optype: moptype_et,

    /// Type of refinement objective
    pub objtype: mobjtype_et,

    /// Controls the debugging output of the program
    pub dbglvl: mdbglvl_et,

    /// The type of coarsening
    pub ctype: mctype_et,

    /// The type of initial partitioning
    pub iptype: miptype_et,

    /// The type of refinement
    pub rtype: mrtype_et,

    /// The # of vertices in the coarsest graph
    pub CoarsenTo: idx_t,

    /// The number of initial partitions to compute
    pub nIparts: idx_t,

    /// Indicates if 2-hop matching will be used
    pub no2hop: idx_t,

    /// Indicates out-of-core execution
    pub ondisk: idx_t,

    /// Indicates if the subdomain connectivity will be minimized
    pub minconn: idx_t,

    /// Indicates if contiguous partitions are required
    pub contig: idx_t,

    /// The number of separators to be found during multiple bisections
    pub nseps: idx_t,

    /// The user-supplied load imbalance factor
    pub ufactor: idx_t,

    /// If the graph will be compressed prior to ordering
    pub compress: idx_t,

    /// If connected components will be ordered separately
    pub ccorder: idx_t,

    /// The seed for the random number generator
    pub seed: idx_t,

    /// The number of different partitionings to compute
    pub ncuts: idx_t,

    /// The number of iterations during each refinement
    pub niter: idx_t,

    /// The user-supplied numflag for the graph
    pub numflag: idx_t,

    /// Indicates if edges will be randomly dropped during coarsening
    pub dropedges: idx_t,

    /// The maximum allowed weight for a vertex
    pub maxvwgt: *mut idx_t,

    /// The number of balancing constraints
    pub ncon: idx_t,

    /// The number of partitions
    pub nparts: idx_t,

    /// .1*(user-supplied prunning factor)
    pub pfactor: real_t,

    /// The per-constraint unbalance factors
    pub ubfactors: *mut real_t,

    /// The target partition weights
    pub tpwgts: *mut real_t,

    ///  The nparts*ncon multiplies for the ith partition and jth constraint for obtaining the balance
    pub pijbm: *mut real_t,

    /// The achieved compression factor
    pub cfactor: real_t,

    /* Various Timers */
    pub TotalTmr: f64,
    pub InitPartTmr: f64,
    pub MatchTmr: f64,
    pub ContractTmr: f64,
    pub CoarsenTmr: f64,
    pub UncoarsenTmr: f64,
    pub RefTmr: f64,
    pub ProjectTmr: f64,
    pub SplitTmr: f64,
    pub Aux1Tmr: f64,
    pub Aux2Tmr: f64,
    pub Aux3Tmr: f64,

    /// Workspace Information:
    ///
    /// The persistent memory core for within function mallocs/frees
    pub mcore: *mut gk_mcore_t,

    /// For use by the k-way refinement routines
    ///
    /// The maximum number of {c,v}nbr_t entries that will ever be allocated
    pub nbrpoolsize_max: usize,

    /// For use by the k-way refinement routines
    ///
    /// The number of {c,v}nbr_t entries that have been allocated
    pub nbrpoolsize: usize,

    /// For use by the k-way refinement routines
    ///
    /// The position of the first free entry in the array
    pub nbrpoolcpos: usize,

    /// For use by the k-way refinement routines
    ///
    /// The number of times the pool was resized
    pub nbrpoolreallocs: usize,

    /// The pool of cnbr_t entries to be used during refinement. The size and current position of the pool is controlled by nnbrs & cnbrs
    pub cnbrpool: *mut cnbr_t,

    /// The pool of vnbr_t entries to be used during refinement. The size and current position of the pool is controlled by nnbrs & cnbrs
    pub vnbrpool: *mut vnbr_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// The maximum allocated number of adjacent domains
    pub maxnads: *mut idx_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// The number of adjacent domains
    pub nads: *mut idx_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// The IDs of the adjacent domains
    pub adids: *mut *mut idx_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// The edge-weight to the adjacent domains
    pub adwgts: *mut *mut idx_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// Auxiliary nparts-size vectors for efficiency
    pub pvec1: *mut idx_t,

    /// Components of the subdomain graph, in sparse format:
    ///
    /// Auxiliary nparts-size vectors for efficiency
    pub pvec2: *mut idx_t,

    /// ondisk related info
    pub pid: libc::pid_t,
}

#[repr(C)]
pub struct graph_t {
    /// The # of vertices and edges in the graph
    pub nvtxs: idx_t,

    /// The # of vertices and edges in the graph
    pub nedges: idx_t,

    /// The # of constrains
    pub ncon: idx_t,

    /// Pointers to the locally stored vertices
    pub xadj: *mut idx_t,

    /// Vertex weights
    pub vwgt: *mut idx_t,

    /// Vertex sizes for min-volume formulation
    pub vsize: *mut idx_t,

    /// Array that stores the adjacency lists of nvtxs
    pub adjncy: *mut idx_t,

    /// Array that stores the weights of the adjacency lists
    pub adjwgt: *mut idx_t,

    /// The sum of the vertex weights in the graph
    pub tvwgt: *mut idx_t,

    /// The inverse of the sum of the vertex weights in the graph
    pub invtvwgt: *mut real_t,

    /// This are to keep track control if the corresponding fields correspond to application or
    /// library memory
    pub free_xadj: std::ffi::c_int,

    /// This are to keep track control if the corresponding fields correspond to application or
    /// library memory
    pub free_vwgt: std::ffi::c_int,

    /// This are to keep track control if the corresponding fields correspond to application or
    /// library memory
    pub free_vsize: std::ffi::c_int,

    /// This are to keep track control if the corresponding fields correspond to application or
    /// library memory
    pub free_adjncy: std::ffi::c_int,

    /// This are to keep track control if the corresponding fields correspond to application or
    /// library memory
    pub free_adjwgt: std::ffi::c_int,

    /// The contraction/coarsening map
    ///  
    /// Gavin: idk the order but it's probably cmap[i] = coarser[i] and so cmap[i] < coarser.nvtxs
    pub cmap: *mut idx_t,

    /// The labels of the vertices for recursive bisection (pmetis/ometis)
    pub label: *mut idx_t,

    /// Partition parameters
    pub mincut: idx_t,

    /// Partition parameters
    pub minvol: idx_t,

    /// Partition parameters
    pub where_: *mut idx_t,

    /// Partition parameters
    ///
    /// Partition weights
    pub pwgts: *mut idx_t,

    /// Partition parameters
    ///
    /// size of boundary for Dal (see appendix)
    pub nbnd: idx_t,

    /// Partition parameters
    ///
    /// lptr of Dal (see appendix)
    pub bndptr: *mut idx_t,

    /// Partition parameters
    ///
    /// lind of Dal (see appendix)
    pub bndind: *mut idx_t,

    /* Bisection refinement parameters */
    pub id: *mut idx_t,
    pub ed: *mut idx_t,

    /// K-way refinement parameter:
    ///
    /// The per-vertex cut-based refinement info
    pub ckrinfo: *mut ckrinfo_t,
    /// K-way refinement parameter:
    ///
    /// The per-vertex volume-based refinement info
    pub vkrinfo: *mut vkrinfo_t,

    /// Node refinement information
    pub nrinfo: *mut nrinfo_t,

    /// various fields for out-of-core processing
    pub gID: std::ffi::c_int,
    /// various fields for out-of-core processing
    pub ondisk: std::ffi::c_int,

    /// keep track of the dropped edgewgt
    pub droppedewgt: idx_t,

    /// the linked-list structure of the sequence of graphs
    pub coarser: *mut graph_t,
    /// the linked-list structure of the sequence of graphs
    pub finer: *mut graph_t,
}

/// The following data structure holds information on degrees for k-way vol-based partition
#[repr(C)]
#[derive(Default, Clone)]
pub struct vkrinfo_t {
    /// The internal degree of a vertex (count of edges)
    pub nid: idx_t,

    /// The total external degree of a vertex (count of edges)
    pub ned: idx_t,

    /// The volume gain of moving that vertex
    pub gv: idx_t,

    /// The number of neighboring subdomains
    pub nnbrs: idx_t,

    /// The index in the vnbr_t array where the nnbrs list of neighbors is stored
    pub inbr: idx_t,
}

/// The following data structure stores holds information on degrees for k-way partition
#[repr(C)]
#[derive(Default)]
pub struct ckrinfo_t {
    /// The internal degree of a vertex (sum of weights)
    pub id: idx_t,

    /// The total external degree of a vertex
    pub ed: idx_t,

    /// The number of neighboring subdomains
    pub nnbrs: idx_t,

    /// The index in the cnbr_t array where the nnbrs list of neighbors is stored
    pub inbr: idx_t,
}

/// The following data structure holds information on degrees for k-way partition
#[repr(C)]
pub struct nrinfo_t {
    pub edegrees: [idx_t; 2],
}

/// This data structure stores volume-based k-way refinement info about an
/// adjacent subdomain for a given vertex.
///
/// Gavin: I believe it stands for volume neighborhood
#[repr(C)]
#[derive(Default, Clone)]
pub struct vnbr_t {
    /// The partition ID
    pub pid: idx_t,

    /// The number of the adjacent edges that are incident on pid
    pub ned: idx_t,

    /// The gain in volume achieved by moving the vertex to pid
    pub gv: idx_t,
}

/// This data structure stores cut-based k-way refinement info about an
/// adjacent subdomain for a given vertex.
#[derive(Clone, Copy)]
#[repr(C)]
pub struct cnbr_t {
    /// The partition ID
    pub pid: idx_t,

    /// The sum of the weights of the adjacent edges
    pub ed: idx_t,
}
