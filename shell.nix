{ pkgs ? import <nixpkgs> {} }:
let
  overrides = (builtins.fromTOML (builtins.readFile ./rust-toolchain.toml));
  libPath = with pkgs; lib.makeLibraryPath [
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
      "-DASSERT2=ON"
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
      "-DASSERT=1"
      "-DASSERT2=1"
    ];
  });
in
  pkgs.mkShell rec {
    name = "metis";

    buildInputs = with pkgs; [
      clang
      llvmPackages_latest.bintools
      rustup
      guile
      new_metis
      cloc
    ];
    RUSTC_VERSION = overrides.toolchain.channel;
    # https://github.com/rust-lang/rust-bindgen#environment-variables
    LIBCLANG_PATH = pkgs.lib.makeLibraryPath [ pkgs.llvmPackages_latest.libclang.lib ];
    shellHook = ''
      export PATH=$PATH:''${CARGO_HOME:-~/.cargo}/bin
      export PATH=$PATH:''${RUSTUP_HOME:-~/.rustup}/toolchains/$RUSTC_VERSION-x86_64-unknown-linux-gnu/bin/
      '';
    # Add precompiled library to rustc search path
    RUSTFLAGS = (builtins.map (a: ''-L ${a}/lib'') [
      # add libraries here (e.g. pkgs.libvmi)
    ]) ++ [
        "-Zbox-noalias=false"
        "-Zmutable-noalias=false"
      ];
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


    hardeningDisable = ["bindnow"];
  }
