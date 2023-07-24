fn main() {
    println!("cargo:rustc-cfg=bindings_module");
    println!("cargo:rustc-link-search=../build");
    println!("cargo:rustc-link-search=../GKlib/build/install/lib");
    println!("cargo:rustc-link-lib=macros");
    println!("cargo:rustc-link-lib=metis");
    println!("cargo:rustc-link-lib=GKlib");
}
