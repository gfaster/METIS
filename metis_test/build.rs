fn main() {
    println!("cargo:rustc-link-search=../libmetis");
    println!("cargo:rustc-link-lib=metis");
}
