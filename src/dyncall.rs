use std::{collections::{HashMap, HashSet}, ffi::{c_void, CStr, CString}, ptr::NonNull, sync::{LazyLock, OnceLock}};

macro_rules! make_cstr {
    ($s:expr) => {{
        const BASE: &str = $s;
        const LEN: usize = BASE.len() + 1;
        const RET_P: [u8; LEN] = const {
            let mut ret: [u8; LEN] = [0; LEN];
            let mut idx = 0;
            loop {
                if idx == LEN - 1 {
                    break;
                }
                ret[idx] = BASE.as_bytes()[idx];
                idx += 1;
            }
            ret
        };
        const { unsafe { std::ffi::CStr::from_bytes_with_nul_unchecked(&RET_P) } }
    }};
}

pub static LIBMETIS: Library = Library::new(
    make_cstr!(env!("LIBMETIS_PORTED"))
);

pub struct Library {
    name: &'static CStr,
    handle: OnceLock<NonNull<c_void>>,
}

impl Library {
    pub const fn new(name: &'static CStr) -> Self {
        Library {
            name,
            handle: OnceLock::new(),
        }
    }

    // pub const fn new(name: &'static CStr) -> Self {
    //     Library {
    //         name,
    //         handle: OnceLock::new(),
    //     }
    // }

    #[inline(never)]
    fn get(&'static self) -> NonNull<c_void> {
        *self.handle.get_or_init(|| {
            // it's important that we use RTLD_LAZY since that works around symbol resolution
            // order. Since the dynamic dispatch within the C version utilizes GNU ifuncs, there's
            // a whole bunch of extra rules.
            let h = unsafe { libc::dlopen(self.name.as_ptr(), libc::RTLD_LAZY) };
            let Some(h) = NonNull::new(h) else {
                let err = unsafe { libc::dlerror() };
                if let Some(err) = NonNull::new(err) {
                    let err = unsafe { CStr::from_ptr(err.as_ptr()) };
                    let err = err.to_string_lossy().to_owned();
                    panic!("Error opening library: {err}");
                } else {
                    panic!("Error opening library: unknown");
                }
            };
            h
        })
    }
}

unsafe impl Send for Library {}
unsafe impl Sync for Library {}

/// exported symbols that are never prefixed
const EXPORTS: &[&str] = &[
];

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
enum Version {
    Rust,
    C
}

const VAR: &str = "METIS_OVERRIDE_SYMS";
static SYM_OVERRIDES: LazyLock<HashMap<CString, Version>> = LazyLock::new(init_overrides);
fn init_overrides() -> HashMap<CString, Version> {
    use std::io::Write;
    let Some(args) = std::env::var_os(VAR) else {
        return HashMap::new()
    };
    let Ok(args) = args.into_string() else {
        let mut out = std::io::stderr();
        writeln!(out, "{VAR} is invalid utf-8").unwrap();
        return HashMap::new()
    };
    let mut ret = HashMap::new();
    for arg in args.split(',') {
        let (sym, spec) = if let Some(split@(sym, _)) = arg.split_once(':') {
            if sym == "c" || sym == "rs" {
                let mut out = std::io::stderr();
                writeln!(out, "Schema: <symbol>:<version> OR <full_symbol>").unwrap();
                continue
            }
            split
        } else if let Some(sym) = arg.strip_prefix("c__") {
            (sym, "c")
        } else if let Some(sym) = arg.strip_prefix("rs__") {
            (sym, "rs")
        } else {
            (arg, "c")
        };
        let lib_pfx = if sym.starts_with("libmetis__") || EXPORTS.contains(&sym) {
            ""
        } else {
            "libmetis__"
        };
        let ver = {
            match spec {
                "c" => Version::C,
                "rs" => Version::Rust,
                _ => {
                    let mut out = std::io::stderr();
                    writeln!(out, "Bad spec: {spec:?}").unwrap();
                    continue
                }
            }
        };
        // always c__ since that's what we lookup with dlsym
        let sym = format!("c__{lib_pfx}{sym}\0");
        let sym = CString::from_vec_with_nul(sym.into_bytes()).unwrap();
        ret.insert(sym, ver);
    }
    ret
}

fn clear_dlerror() {
    unsafe { libc::dlerror() };
}

fn get_dlerror() -> Option<String> {
    let p = NonNull::new(unsafe { libc::dlerror() })?;
    // Note the on glibc and more recent musl, dlerror is thread safe:
    // https://wiki.musl-libc.org/functional-differences-from-glibc.html
    let cstr = unsafe { CStr::from_ptr(p.as_ptr()) };
    let s = cstr.to_string_lossy().to_string();
    Some(s)
}

pub struct ICall {
    func: OnceLock<NonNull<c_void>>,
    sym_name: &'static CStr,
    library: Option<&'static Library>,
    rs_ver: NonNull<c_void>,
    c_ver: Option<NonNull<c_void>>,
}

unsafe impl Send for ICall {}
unsafe impl Sync for ICall {}

impl ICall {
    #[allow(dead_code)]
    pub const unsafe fn new(c_sym: &'static CStr, library: &'static Library, rs_ver: NonNull<c_void>) -> Self {
        ICall {
            func: OnceLock::new(),
            sym_name: c_sym,
            library: Some(library),
            rs_ver,
            c_ver: None,
        }
    }

    #[allow(dead_code)]
    pub const unsafe fn new_referenced(c_sym: &'static CStr, c_ver: NonNull<c_void>, rs_ver: NonNull<c_void>) -> Self {
        ICall {
            func: OnceLock::new(),
            sym_name: c_sym,
            library: None,
            rs_ver,
            c_ver: Some(c_ver),
        }
    }

    pub fn get(&'static self) -> NonNull<c_void> {
        // let overrides = &*SYM_OVERRIDES;
        // println!("{overrides:?}");
        // panic!("");
        *self.func.get_or_init(|| {
            let ver = if let Some(&ver) = SYM_OVERRIDES.get(self.sym_name) {
                ver
            } else {
                Version::Rust
            };
            match ver {
                Version::Rust => self.rs_ver,
                Version::C => {
                    if let Some(c_ver) = self.c_ver {
                        c_ver
                    } else {
                        let lib = self.library.expect("library handle object exists").get();
                        clear_dlerror();
                        let sym = unsafe { libc::dlsym(lib.as_ptr(), self.sym_name.as_ptr()) };
                        let Some(sym) = NonNull::new(sym) else {
                            // probably an error state -- but definitely unusable
                            let msg = get_dlerror().unwrap_or_else(|| "Could not get valid handle, but no error".into());
                            panic!("could not resolve `{}`: {msg}", self.sym_name.to_string_lossy());
                        };
                        sym
                    }
                },
            }
        })
    }
}
