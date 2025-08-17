{ pkgs ? import <nixpkgs> {} }:
let
  aux-normalization = pkgs.callPackage ./normalization/aux-normalization.nix { };
  metis-normalized = pkgs.callPackage ./metis-normalized.nix { inherit aux-normalization; };

  overrides = (builtins.fromTOML (builtins.readFile ./rust-toolchain.toml));
  libPath = with pkgs; lib.makeLibraryPath [
    aux-normalization
    # load external libraries that you need in your rust project here
  ];

  gklib = pkgs.stdenv.mkDerivation rec {
    pname = "gklib";
    version = "dev";

    src = pkgs.fetchFromGitHub {
      owner = "KarypisLab";
      repo = "GKlib";
      rev = "8bd6bad750b2b0d90800c632cf18e8ee93ad72d7";
      hash = "sha256-tunepMLaRDR5FQVL/9S7/w6e1j+f2+pg01H/0/z/ZCI=";
    };

    patches = [
      ./gk_fifo.patch
    ];

    cmakeFlags = [
      "-DCMAKE_C_FLAGS_RELEASE=-O2"
      "-DASSERT=ON"
    ];

    nativeBuildInputs = with pkgs; [ cmake ];
  };

  new_metis = pkgs.metis.overrideAttrs (finalAttrs: previousAttrs: rec {
    src = pkgs.fetchFromGitHub {
      owner = "KarypisLab";
      repo = "METIS";
      rev = "v5.2.1";
      hash = "sha256-eddLR6DvZ+2LeR0DkknN6zzRvnW+hLN2qeI+ETUPcac=";
    };

    preConfigure = ''
      mkdir -p build/xinclude
      echo "#define IDXTYPEWIDTH 32" > build/xinclude/metis.h
      echo "#define REALTYPEWIDTH 32" >> build/xinclude/metis.h
      cat include/metis.h >> build/xinclude/metis.h
      cp include/CMakeLists.txt build/xinclude
    '';

    cmakeFlags = [
      "-DCMAKE_C_FLAGS_RELEASE=-O2"
      "-DCMAKE_C_FLAGS=-w"
      "-DGKLIB_PATH=${gklib}"
      "-DASSERT=ON"
    ];
  });

in
  pkgs.mkShell rec {
    name = "metis";

    METIS_NORM = metis-normalized;
    METIS_NORM_SRC = metis-normalized.src;

    METIS_NORM_FUNCTIONS = pkgs.lib.makeLibraryPath [ aux-normalization ];

    buildInputs = with pkgs; [
      clang
      llvmPackages_latest.bintools
      rustup
      new_metis

      cloc
      valgrind
      guile
      ocamlPackages.magic-trace linuxKernel.packages.linux_5_10.perf
    ];
    RUSTC_VERSION = overrides.toolchain.channel;
    # https://github.com/rust-lang/rust-bindgen#environment-variables
    LIBCLANG_PATH = pkgs.lib.makeLibraryPath [ pkgs.llvmPackages_latest.libclang.lib ];

    # See: https://sourceware.org/gdb/current/onlinedocs/gdb.html/File-Options.html
    # See: https://nixos.wiki/wiki/Debug_Symbols
    shellHook = ''
      export PATH=$PATH:''${CARGO_HOME:-~/.cargo}/bin
      export PATH=$PATH:''${RUSTUP_HOME:-~/.rustup}/toolchains/$RUSTC_VERSION-x86_64-unknown-linux-gnu/bin/

      alias norm-gdb="gdb -d ${metis-normalized.src}/programs -d \
        ${metis-normalized.src}/libmetis -d ${metis-normalized.src}/include -x gdb/options_pp.scm"
      '';

    LD_LIBRARY_PATH = libPath;
    # Add glibc, clang, glib, and other headers to bindgen search path
    BINDGEN_EXTRA_CLANG_ARGS = 
      # Includes normal include path
      (builtins.map (a: ''-I"${a}/include"'') [
        # add dev libraries here (e.g. pkgs.libvmi.dev)
        pkgs.glibc.dev
      ])
      # Includes with special directory paths
      ++ [
        ''-I"${pkgs.llvmPackages_latest.libclang.lib}/lib/clang/${pkgs.llvmPackages_latest.libclang.version}/include"''
        ''-I"${pkgs.glib.dev}/include/glib-2.0"''
        ''-I${pkgs.glib.out}/lib/glib-2.0/include/''
      ];

    # needed for dyncall system
    hardeningDisable = ["bindnow"];
  }
