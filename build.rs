use std::{fs::read_dir, os::unix::prelude::OsStrExt};

use cc::Build;

fn main() {
    let mut files: Vec<_> = read_dir("src")
        .expect("src/ files exist")
        .filter_map(|f| f.ok())
        .filter(|f| f.file_name().as_bytes().ends_with(b".c"))
        .map(|f| f.path())
        .collect();

    // dbg!(std::env::vars().filter(|(key, _)| key.contains("CARGO")).collect::<Vec<_>>());
    let do_dual_link: i32 = std::env::var("CARGO_FEATURE_DUAL_LINK")
        .unwrap_or("0".to_string())
        .parse()
        .expect("cargo feature is either 1 or 0");

    let ported_files: Vec<_> = read_dir("src/ported")
        .expect("src/ported files exist")
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

    println!("cargo:rerun-if-changed=src/ported");
    println!("cargo:rerun-if-changed=build.rs");

    if do_dual_link != 0 {
        files.extend(dbg!(ported_files));
        // panic!()
    }

    // panic!();

    // dbg!(&files);

    Build::new()
        .files(files)
        .compiler("gcc")
        .include("GKlib/build/install/include")
        .include("include")
        .include("src")
        .define("IDXTYPEWIDTH", "32")
        .define("REALTYPEWIDTH", "32")
        .define("DMALLOC", "")
        // .define("ASSERT", "1")
        // .define("ASSERT2", "1")
        .warnings(false)
        // .pic(false)
        // .link_lib_modifier("-bundle")
        .compile("metis");

    println!("cargo:rustc-link-search=GKlib/build/install/lib");
    println!("cargo:rustc-link-lib=GKlib");
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=GKlib");
}
