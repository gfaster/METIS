use std::{collections::{HashMap, HashSet}, ffi::{c_void, CStr, CString}, ptr::NonNull, sync::{LazyLock, OnceLock}};

pub static LIBMETIS: Library = Library::new(c"");

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

    #[inline(never)]
    fn get(&'static self) -> NonNull<c_void> {
        *self.handle.get_or_init(|| {
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
        let Some((sym, spec)) = arg.split_once(':') else {
            let mut out = std::io::stderr();
            writeln!(out, "Not split").unwrap();
            continue
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
        let sym = format!("{lib_pfx}{sym}\0");
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
    library: &'static Library,
    rs_ver: NonNull<c_void>,
}

unsafe impl Send for ICall {}
unsafe impl Sync for ICall {}

impl ICall {
    pub const unsafe fn new(c_sym: &'static CStr, library: &'static Library, rs_ver: NonNull<c_void>) -> Self {
        ICall {
            func: OnceLock::new(),
            sym_name: c_sym,
            library,
            rs_ver,
        }
    }

    pub fn get(&'static self) -> NonNull<c_void> {
        *self.func.get_or_init(|| {
            let ver = if let Some(&ver) = SYM_OVERRIDES.get(self.sym_name) {
                ver
            } else {
                Version::C
            };
            match ver {
                Version::Rust => self.rs_ver,
                Version::C => {
                    let lib = self.library.get();
                    clear_dlerror();
                    let sym = unsafe { libc::dlsym(lib.as_ptr(), self.sym_name.as_ptr()) };
                    let Some(sym) = NonNull::new(sym) else {
                        // probably an error state -- but definitely unusable
                        let msg = get_dlerror().unwrap_or_else(|| "Could not get valid handle, but no error".into());
                        panic!("could not resolve `{}`: {msg}", self.sym_name.to_string_lossy());
                    };
                    sym
                },
            }
        })
    }
}
