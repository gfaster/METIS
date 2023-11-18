use std::{fs::read_dir, os::unix::prelude::OsStrExt, process::Command};

use cc::Build;

fn main() {
    let do_no_rs = std::env::var("CARGO_FEATURE_NO_RS")
        .unwrap_or("0".to_string())
        .parse::<i32>()
        .expect("cargo feature is either 1 or 0")
        != 0;

    let mut files: Vec<_> = read_dir("src")
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
    let do_dual_link = std::env::var("CARGO_FEATURE_DUAL_LINK")
        .unwrap_or("0".to_string())
        .parse::<i32>()
        .expect("cargo feature is either 1 or 0")
        != 0;

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

    for file in &gklib_files {
        println!(
            "cargo:rerun-if-changed={}",
            file.to_str().expect("files are valid utf-8")
        );
    }

    println!("cargo:rerun-if-changed=src/ported");
    println!("cargo:rerun-if-changed=build.rs");

    if do_dual_link || do_no_rs {
        files.extend(dbg!(ported_files));
        // panic!()
    }

    // panic!();
    Command::new("pwd").spawn().unwrap().wait().unwrap();

    // dbg!(&files);
    Build::new()
        .files(gklib_files)
        .compiler("gcc")
        .include("GKlib")
        .define("LINUX", "")
        .define("NDEBUG2", "")
        .define("_FILE_OFFSET_BITS", "64")
        .flag("-fno-strict-aliasing")
        .flag("-std=c99")
        .warnings(false)
        .compile("GKlib");

    Build::new()
        .files(files)
        .compiler("gcc")
        .include("GKlib/")
        .include("include")
        .include("src")
        .define("IDXTYPEWIDTH", "32")
        .define("REALTYPEWIDTH", "32")
        .define("DMALLOC", "")
        .define("NDEBUG2", "")
        .flag("-fno-strict-aliasing")
        .warnings(false)
        // .pic(false)
        // .link_lib_modifier("-bundle")
        .compile("metis");

    // println!("cargo:rustc-link-search=GKlib/build/install/lib");
    // println!("cargo:rustc-link-lib=GKlib");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=GKlib");
}
