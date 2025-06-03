//! Allows changing which functions are used (C or Rust) via environment variable
//!
//! Set `METIS_OVERRIDE_SYMS` to do so. See [`translation.md`](../translation.md) for more info.

use crate::util::make_cstr;
use std::{
    borrow::Cow,
    cell::{Cell, RefCell},
    collections::{HashMap, HashSet},
    ffi::{c_void, CStr, CString},
    ptr::NonNull,
    sync::{atomic::AtomicPtr, LazyLock, OnceLock},
    thread::LocalKey,
};

pub static LIBMETIS: Library = Library::new(make_cstr!(env!("LIBMETIS_PORTED")));

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
const EXPORTS: &[&str] = &[];

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
enum Version {
    Rust,
    C,
}

const VAR: &str = "METIS_OVERRIDE_SYMS";
static SYM_OVERRIDES: LazyLock<Overrides> = LazyLock::new(Overrides::init_overrides);
#[derive(Default)]
struct Overrides {
    globs: Vec<(Glob<'static>, Version)>,
    exact: HashMap<Box<[u8]>, Version>,
}

impl Overrides {
    fn get(&self, name: impl AsRef<[u8]>) -> Version {
        let name = name.as_ref();
        let name = name.strip_suffix(&[0u8]).unwrap_or(name);
        if let Some(&exact_ver) = self.exact.get(name) {
            #[cfg(test)]
            ICall::update_last_resolve();

            return exact_ver;
        }
        let short_name = name.strip_prefix(b"c__").unwrap_or(name);
        let short_name = short_name.strip_prefix(b"rs__").unwrap_or(short_name);
        let short_name = short_name.strip_prefix(b"libmetis__").unwrap_or(short_name);
        if let Some(&(_, glob_ver)) = self
            .globs
            .iter()
            .rev()
            .find(|(glob, _)| glob.matches(short_name))
        {
            #[cfg(test)]
            ICall::update_last_resolve();

            return glob_ver;
        }
        Version::Rust
    }

    fn init_with(args: &str) -> Self {
        use std::io::Write;
        let mut ret = Self::default();
        for arg in args.split(',') {
            if arg.contains('*') {
                // this is a glob!
                let (glob, spec) = if let Some(split) = arg.split_once(':') {
                    split
                } else {
                    (arg, "c")
                };
                let ver = {
                    match spec {
                        "c" => Version::C,
                        "rs" => Version::Rust,
                        _ => {
                            let mut out = std::io::stderr();
                            writeln!(out, "Bad spec: {spec:?}").unwrap();
                            continue;
                        }
                    }
                };
                ret.globs.push((Glob::new_owned(glob), ver));
                continue;
            }
            let (sym, spec) = if let Some(split @ (sym, _)) = arg.split_once(':') {
                if sym == "c" || sym == "rs" {
                    let mut out = std::io::stderr();
                    writeln!(out, "Schema: <symbol>:<version> OR <full_symbol>").unwrap();
                    continue;
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
                        continue;
                    }
                }
            };
            // always c__ since that's what we lookup with dlsym
            let sym = format!("c__{lib_pfx}{sym}").into_bytes();
            ret.exact.insert(sym.into(), ver);
        }
        ret
    }

    fn init_overrides() -> Self {
        use std::io::Write;
        let Some(args) = std::env::var_os(VAR) else {
            return Overrides::default();
        };
        let Ok(args) = args.into_string() else {
            let mut out = std::io::stderr();
            writeln!(out, "{VAR} is invalid utf-8").unwrap();
            return Overrides::default();
        };
        Self::init_with(&args)
    }
}

thread_local! {
    static LOCAL_OVERRIDES: RefCell<Overrides> = RefCell::new(Overrides::default());
    // start at 1 so we can allow a dangling pointer and use wrapping-panic for the safety
    static LOCAL_EPOCH: Cell<u64> = const { Cell::new(1) };
    static LOCAL_EPOCH_LAST_RESOLVE: Cell<u64> = const { Cell::new(0) };
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

pub struct OverrideKey {
    _no_construct: (),
}

#[cfg_attr(not(test), expect(dead_code))]
pub fn set_local_overrides(overrides: &str) -> OverrideKey {
    if cfg!(test) {
        println!("\nsetting overrides: `{overrides}`")
    }
    LOCAL_OVERRIDES.with(|l| *l.borrow_mut() = Overrides::init_with(overrides));
    let Some(x) = LOCAL_EPOCH.get().checked_add(1) else {
        panic!("Local epoch overflowed")
    };
    LOCAL_EPOCH.set(x);
    OverrideKey { _no_construct: () }
}

impl Drop for OverrideKey {
    fn drop(&mut self) {
        LOCAL_OVERRIDES.with(|l| *l.borrow_mut() = Overrides::default());
    }
}

pub struct ICall {
    func: OnceLock<NonNull<c_void>>,
    sym_name: &'static CStr,
    library: Option<&'static Library>,
    rs_ver: NonNull<c_void>,
    c_ver: OnceLock<NonNull<c_void>>,
}

unsafe impl Send for ICall {}
unsafe impl Sync for ICall {}

impl ICall {
    #[allow(dead_code)]
    pub const unsafe fn new(
        c_sym: &'static CStr,
        library: &'static Library,
        rs_ver: NonNull<c_void>,
    ) -> Self {
        ICall {
            func: OnceLock::new(),
            sym_name: c_sym,
            library: Some(library),
            rs_ver,
            c_ver: OnceLock::new(),
        }
    }

    #[cfg_attr(not(test), expect(dead_code))]
    fn update_last_resolve() {
        LOCAL_EPOCH_LAST_RESOLVE.set(Self::epoch());
    }

    pub fn epoch() -> u64 {
        LOCAL_EPOCH.get()
    }

    #[cfg_attr(not(test), expect(dead_code))]
    fn last_resolve_epoch() -> u64 {
        LOCAL_EPOCH_LAST_RESOLVE.get()
    }

    #[cfg_attr(not(test), expect(dead_code))]
    pub fn get_dynamic(&'static self) -> NonNull<c_void> {
        let ver = LOCAL_OVERRIDES.with(|ov| ov.borrow().get(self.sym_name.to_bytes()));
        let ret = self.get_ver(ver);
        // eprintln!("resolved {ver:?} of {:?}", self.sym_name);
        ret
    }

    fn get_c_func(&'static self) -> NonNull<c_void> {
        *self.c_ver.get_or_init(|| {
            let lib = self.library.expect("library handle object exists").get();
            clear_dlerror();
            let sym = unsafe { libc::dlsym(lib.as_ptr(), self.sym_name.as_ptr()) };
            let Some(sym) = NonNull::new(sym) else {
                // probably an error state -- but definitely unusable
                let msg = get_dlerror()
                    .unwrap_or_else(|| "Could not get valid handle, but no error".into());
                panic!(
                    "could not resolve `{}`: {msg}",
                    self.sym_name.to_string_lossy()
                );
            };
            sym
        })
    }

    fn get_ver(&'static self, ver: Version) -> NonNull<c_void> {
        match ver {
            Version::Rust => self.rs_ver,
            Version::C => self.get_c_func(),
        }
    }

    #[cfg_attr(test, expect(dead_code))]
    pub fn get(&'static self) -> NonNull<c_void> {
        // let overrides = &*SYM_OVERRIDES;
        // println!("{overrides:?}");
        // panic!("");
        *self.func.get_or_init(|| {
            let ver = SYM_OVERRIDES.get(self.sym_name.to_bytes());
            self.get_ver(ver)
        })
    }
}

/// Helper for grouping functions -- very dumb and can get very slow if there are too many `*`
pub struct Glob<'a> {
    template: Cow<'a, [u8]>,
}

impl Glob<'static> {
    #[allow(dead_code)]
    pub fn new_owned(g: impl AsRef<[u8]>) -> Self {
        Self {
            template: Cow::Owned(g.as_ref().to_owned()),
        }
    }
}

impl<'a> Glob<'a> {
    #[allow(dead_code)]
    pub const fn new_str(b: &'a str) -> Self {
        Self {
            template: Cow::Borrowed(b.as_bytes()),
        }
    }

    #[allow(dead_code)]
    pub const fn new_bytes(b: &'a [u8]) -> Self {
        Self {
            template: Cow::Borrowed(b),
        }
    }

    // TODO: optimize me!
    pub fn matches(&self, s: impl AsRef<[u8]>) -> bool {
        fn slices(s: &[u8]) -> impl Iterator<Item = &[u8]> {
            (0..s.len()).map(|i| &s[i..])
        }
        fn subslices<'a>(haystack: &'a [u8], needle: &'a [u8]) -> impl Iterator<Item = &'a [u8]> {
            slices(haystack).filter_map(|slice| slice.strip_prefix(needle))
        }
        fn initial(mut g: &[u8], mut s: &[u8]) -> bool {
            let Some(star_idx) = g.iter().position(|&c| c == b'*') else {
                return g == s;
            };
            if &g[..star_idx] != &s[..star_idx] {
                return false;
            }
            g = &g[star_idx + 1..];
            s = &s[star_idx..];
            inner(g, s)
        }
        /// assumes g starts with an implicit `*`
        fn inner(mut g: &[u8], s: &[u8]) -> bool {
            // eprintln!("called with => g: {:?}, s: {:?}", std::str::from_utf8(g), std::str::from_utf8(s));
            let leading_stars = g.iter().take_while(|&&c| c == b'*').count();
            g = &g[leading_stars..];
            if g.is_empty() {
                // eprintln!("empty glob => {:?}", std::str::from_utf8(s));
                return true;
            }
            if let Some(lit_len) = g.iter().position(|&c| c == b'*') {
                debug_assert!(lit_len >= 1);
                let lit = &g[..lit_len];
                g = &g[lit_len..];
                for s in subslices(s, lit) {
                    if inner(g, s) {
                        return true;
                    }
                }
                false
            } else {
                s.ends_with(g)
            }
        }
        initial(&*self.template, s.as_ref())
    }
}

#[cfg(test)]
fn floor_char_boundary(s: &str, idx: usize) -> usize {
    if idx >= s.len() {
        return s.len()
    }
    idx - (idx.saturating_sub(3)..=idx).rev().position(|i| s.is_char_boundary(i)).unwrap()
}

#[cfg(test)]
#[track_caller]
fn ab_assert_eq<T>(l: T, r: T) 
    where T: Eq + std::fmt::Debug 
{
    use std::fmt::Write;
    if l != r {
        let mut buf1 = String::with_capacity(1024);
        write!(buf1, "{l:?}").unwrap();
        let i1 = floor_char_boundary(&buf1, 1024);
        let s1 = &buf1[..i1];
        let dot1 = if buf1.len() > i1 {
            "..."
        } else {
            ""
        };


        let mut buf2 = String::with_capacity(1024);
        write!(buf2, "{r:?}").unwrap();
        let i2 = floor_char_boundary(&buf2, 1024);
        let s2 = &buf2[..i2];
        let dot2 = if buf2.len() > i2 {
            "..."
        } else {
            ""
        };

        panic!("A/B test not equal!\n  \
            left: {s1}{dot1}\n  \
            right: {s2}{dot2}\n
            ")
    }
}

/// runs the provided function twice. Once with the first overrides and once with the second overrides
///
/// Will panic if the provided override isn't used
#[cfg(test)]
pub(crate) fn ab_test_multi<F, T>(overrides: (&str, &str), mut f: F) -> (T, T)
where
    F: FnMut() -> T,
{
    let seed = 4321;
    fastrand::seed(seed);
    let key = crate::dyncall::set_local_overrides(overrides.0);
    let no_override = f();
    drop(key);
    let key = crate::dyncall::set_local_overrides(overrides.1);
    fastrand::seed(seed);
    let with_overrides = f();
    drop(key);
    let last_resolve_epoch = ICall::last_resolve_epoch();
    let current_epoch = ICall::epoch();
    assert!(
        last_resolve_epoch == current_epoch,
        "No dynamic overrides actually occurred in AB test!"
    );
    (no_override, with_overrides)
}

/// runs the provided function twice. Once with no overrides and once with the overrides specified
/// by the overrides argument
///
/// Will panic if the provided override isn't used
#[cfg(test)]
#[expect(unused)]
pub(crate) fn ab_test_multi_eq<F, T>(overrides: (&str, &str), f: F)
where
    F: FnMut() -> T,
    T: Eq + std::fmt::Debug,
{
    let (no_override, with_overrides) = ab_test_multi(overrides, f);
    ab_assert_eq(no_override, with_overrides)
}

/// runs the provided function twice. Once with no overrides and once with the overrides specified
/// by the overrides argument
///
/// Will panic if the provided override isn't used
#[cfg(test)]
pub(crate) fn ab_test<F, T>(overrides: &str, f: F) -> (T, T)
where
    F: FnMut() -> T,
{
    ab_test_multi(("", overrides), f)
}

/// runs the provided function three times. the second time with all C functions and the third with all C
/// functions, with overrides added, returning the second and third runs. The first run is to
/// ensure the override is used
#[cfg(test)]
pub(crate) fn ab_test_single<F, T>(overrides: &str, mut f: F) -> (T, T)
where
    F: FnMut() -> T,
{
    let exc_overrides = format!("*,{overrides}");
    let ret = ab_test_multi(("*", &exc_overrides), &mut f);

    {
        // make sure the override is used - comes after to assist in debugging crashes
        fastrand::seed(0);
        let key = crate::dyncall::set_local_overrides(overrides);
        let _ = f();
        drop(key);
        let last_resolve_epoch = ICall::last_resolve_epoch();
        let current_epoch = ICall::epoch();
        assert!(
            last_resolve_epoch == current_epoch,
            "No dynamic overrides actually occurred in AB test!"
        );
    }
    ret
}

/// runs the provided function three times. the second time with all C functions and the third with all C
/// functions, with overrides added. Asserts the results of the second and third are equal
///
/// Will panic if the provided override isn't used
#[cfg(test)]
pub(crate) fn ab_test_single_eq<F, T>(overrides: &str, f: F)
where
    F: FnMut() -> T,
    T: Eq + std::fmt::Debug,
{
    let (no_override, with_overrides) = ab_test_single(overrides, f);
    ab_assert_eq(no_override, with_overrides)
}

/// runs the provided function twice. Once with no overrides and once with the overrides specified
/// by the overrides argument
///
/// Will panic if the provided override isn't used
#[cfg(test)]
#[track_caller]
pub(crate) fn ab_test_eq<F, T>(overrides: &str, f: F)
where
    F: FnMut() -> T,
    T: Eq + std::fmt::Debug,
{
    let (no_override, with_overrides) = ab_test(overrides, f);
    ab_assert_eq(no_override, with_overrides)
}

#[cfg(test)]
#[metis_func]
pub unsafe extern "C" fn TestVerify() -> u32 {
    0xdeadbeaf
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glob_matches() {
        #[track_caller]
        fn case(glob: &str, haystack: &str) {
            eprintln!("begin with ===> g: {glob:?}, s: {haystack:?}");
            let g = Glob::new_bytes(glob.as_bytes());
            assert!(
                g.matches(haystack),
                "glob {glob:?} did not match {haystack:?}"
            )
        }
        case("abc", "abc");
        case("a", "a");
        case("", "");
        case("*", "abc");
        case("*", "a");
        case("*", "");
        case("a*", "a");
        case("a*", "abc");
        case("*a*", "abc");
        case("*A*", "  A  ");
        case("*A", "  A");
        case("*A", "  AA");
        case("*A", "  CBA");
        case("*ABC", "abcABC");
        case("S*MID*E", "S---MID---E");
        case("S**MID**E", "S---MID---E");
        case("S*1*2*E", "S---1--2---E");
        case("S*1*2*E", "S---12---E");
        case("S*12*3*E", "S---12--3---E");
    }

    #[test]
    fn glob_matches_not() {
        #[track_caller]
        fn case(glob: &str, haystack: &str) {
            let g = Glob::new_bytes(glob.as_bytes());
            assert!(!g.matches(haystack), "glob {glob:?} matched {haystack:?}")
        }
        case("", "abc");
        case("a", "ab");
        case("a", "ba");
        case("a*", "ba");
        case("a*", "bac");
        case("*A", "--A-");
        case("*A", "--AA-");
        case("*A", "A-");
        case("S*1*2*E", "S---13---E");
    }

    #[test]
    fn runtime_version_switching() {
        const RUST_EXPECTED: u32 = 0xdeadbeaf;
        const C_EXPECTED: u32 = 0xabad1dea;

        let (rs_ver, c_ver) = ab_test("TestVerify", || unsafe { TestVerify() });

        assert_eq!(
            (rs_ver, c_ver),
            (0xdeadbeaf, 0xabad1dea),
            "Runtime version switching failed:\n  \
            rust version:\n    \
            expected: {RUST_EXPECTED:#X}\n    \
            actual: {rs_ver:#X}\n  \
            c version:\n    \
            expected: {C_EXPECTED:#X}\n    \
            actual: {c_ver:#X}\n  \
            "
        )
    }

    #[test]
    fn runtime_version_switching_from_c() {
        const RUST_EXPECTED: u32 = 0xdeadbeaf;
        const C_EXPECTED: u32 = 0xabad1dea;

        extern "C" {
            fn call_test_verify_from_c() -> u32;
        }

        let (rs_ver, c_ver) = ab_test("TestVerify", || unsafe { call_test_verify_from_c() });

        assert_eq!(
            (rs_ver, c_ver),
            (0xdeadbeaf, 0xabad1dea),
            "Runtime version switching failed:\n  \
            rust version:\n    \
            expected: {RUST_EXPECTED:#X}\n    \
            actual: {rs_ver:#X}\n  \
            c version:\n    \
            expected: {C_EXPECTED:#X}\n    \
            actual: {c_ver:#X}\n  \
            "
        )
    }

    #[test]
    #[should_panic = "No dynamic overrides"]
    fn ab_test_fail_when_bad_symbol() {
        // oh no, we misspelled it
        let _ = ab_test("TetVerify", || unsafe { TestVerify() });
    }
}
