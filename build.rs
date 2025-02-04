use std::{env, fs::read_dir, os::unix::prelude::OsStrExt, path::PathBuf, process::Command};

use cc::Build;

trait PassThrough {
    fn pass(&mut self, f: impl FnOnce(&mut Self)) -> &mut Self {
        f(self);
        self
    }
}
impl PassThrough for Build {}

fn main() {
    // let do_no_rs = std::env::var("CARGO_FEATURE_NO_RS")
    //     .unwrap_or("0".to_string())
    //     .parse::<i32>()
    //     .expect("cargo feature is either 1 or 0")
    //     != 0;

    let files: Vec<_> = read_dir("src")
        .expect("src/ files exist")
        .filter_map(|f| f.ok())
        .filter(|f| f.file_name().as_bytes().ends_with(b".c"))
        .map(|f| f.path())
        .collect();

    for file in &files {
        println!(
            "cargo:rerun-if-changed={}",
            file.to_str().expect("files are valid utf-8")
        );
    }

    // dbg!(std::env::vars().filter(|(key, _)| key.contains("CARGO")).collect::<Vec<_>>());
    // let do_dual_link = std::env::var("CARGO_FEATURE_DUAL_LINK")
    //     .unwrap_or("0".to_string())
    //     .parse::<i32>()
    //     .expect("cargo feature is either 1 or 0")
    //     != 0;

    let ported_files: Vec<_> = read_dir("src/ported")
        .expect("src/ported files exist")
        .filter_map(|f| f.ok())
        .filter(|f| f.file_name().as_bytes().ends_with(b".c"))
        .map(|f| f.path())
        .collect();

    let mut headers: Vec<_> = read_dir("src")
        .expect("src/ files exist")
        .filter_map(|f| f.ok())
        .filter(|f| f.file_name().as_bytes().ends_with(b".h"))
        .map(|f| f.path())
        .collect();
    headers.extend(
        read_dir("GKlib")
            .expect("GKlib/ files exist")
            .filter_map(|f| f.ok())
            .filter(|f| f.file_name().as_bytes().ends_with(b".h"))
            .map(|f| f.path()),
    );

    for file in &headers {
        println!(
            "cargo:rerun-if-changed={}",
            file.to_str().expect("files are valid utf-8")
        );
    }

    let gklib_files: Vec<_> = read_dir("GKlib")
        .expect("GKlib/ files exist")
        .filter_map(|f| f.ok())
        .filter(|f| f.file_name().as_bytes().ends_with(b".c"))
        .map(|f| f.path())
        .collect();

    assert!(!gklib_files.is_empty(), "found no gklib files");

    for file in &gklib_files {
        println!(
            "cargo:rerun-if-changed={}",
            file.to_str().expect("files are valid utf-8")
        );
    }

    println!("cargo:rerun-if-changed=src/ported");
    println!("cargo:rerun-if-changed=build.rs");

    // files.extend(ported_files.clone());

    // panic!();
    Command::new("pwd").spawn().unwrap().wait().unwrap();

    let ndebug = !cfg!(debug_assertions);

    // dbg!(&files);
    Build::new()
        .files(gklib_files)
        .compiler("gcc")
        .include("GKlib")
        .define("LINUX", "")
        .define("NDEBUG2", "")
        .pass(|b| {
            if ndebug {
                b.define("NDEBUG", "");
            }
        })
        .define("_FILE_OFFSET_BITS", "64")
        .flag("-fno-strict-aliasing")
        .flag("-std=c99")
        .warnings(false)
        .compile("GKlib");

    let mut metis = Build::new();
    metis.files(files)
        .compiler("gcc")
        .include("GKlib/")
        .include("include")
        .include("src")
        .define("IDXTYPEWIDTH", "32")
        .define("REALTYPEWIDTH", "32")
        .define("DMALLOC", "")
        .define("NDEBUG2", "")
        .pass(|b| {
            if ndebug {
                b.define("NDEBUG", "");
            }
        })
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .compile("metis_unported")
        // .pic(false)
        // .link_lib_modifier("-bundle")
        // .compile("metis")
    ;

    make_so(&ported_files);

    // println!("cargo:rustc-link-search=GKlib/build/install/lib");
    println!("cargo:rustc-link-lib=GKlib");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=GKlib");
    println!("cargo:rustc-link-arg=-Wl,--export-dynamic");
}

fn make_so(source: &[PathBuf]) {
    let ndebug = !cfg!(debug_assertions);
    let objs = Build::new().files(source)
        .compiler("gcc")
        .include("GKlib/")
        .include("include")
        .include("src")
        .define("IDXTYPEWIDTH", "32")
        .define("REALTYPEWIDTH", "32")
        .define("DMALLOC", "")
        .define("NDEBUG2", "")
        .pass(|b| {
            if ndebug {
                b.define("NDEBUG", "");
            }
        })
        .flag("-fno-strict-aliasing")
        .warnings(false)
        .compile_intermediates()
        // .pic(false)
        // .link_lib_modifier("-bundle")
        // .compile("metis")
    ;

    let out = env::var("OUT_DIR").unwrap();
    let outfile = format!("{out}/libmetis_ported.so");
    println!("cargo::rustc-env=LIBMETIS_PORTED={outfile}");
    let mut cmd = Command::new("ld")
        .arg("--shared")
        .arg("-o")
        .arg(outfile)
        .args(objs)
        .arg("-L")
        .arg(out)
        .arg("-lc")
        .arg("-lm")
        .arg("-lmetis_unported")
        .arg("-lGKlib")
        .spawn()
        .unwrap();
    let status = cmd.wait().unwrap();
    assert!(status.success())
}
