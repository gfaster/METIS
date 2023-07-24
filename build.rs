use std::{os::unix::prelude::OsStrExt, fs::read_dir};

use cc::Build;

fn main() {
    let files: Vec<_> = read_dir("src").expect("src/ files exist").filter_map(|f| f.ok()).filter(|f| f.file_name().as_bytes().ends_with(b".c")).map(|f| f.path()).collect();

    dbg!(&files);

    Build::new()
        .files(files)
        .compiler("gcc")
        .include("GKlib/build/install/include")
        .include("include")
        .include("src")
        .define("IDXTYPEWIDTH", "32")
        .define("REALTYPEWIDTH", "32")
        // .define("ASSERT", "1")
        // .define("ASSERT2", "1")
        .warnings(false)
        // .pic(false)
        // .link_lib_modifier("-bundle")
        .compile("metis");


    println!("cargo:rustc-link-search=GKlib/build/install/lib");
    println!("cargo:rustc-link-lib=GKlib");

}
