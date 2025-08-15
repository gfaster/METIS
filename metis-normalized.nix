{ pkgs ? import <nixpkgs> {} }:
let
  aux-impls = pkgs.stdenv.mkDerivation rec {
    pname = "functions";
    version = "dev";

    src = ./normalization;

    nativeBuildInputs = with pkgs; [ rustc ];

    buildPhase = ''
      runHook preBuild
      rustc functions.rs
      runHook postBuild
    '';

    installPhase = ''
      runHook preInstall
      mkdir -p $out/lib
      mv libfunctions.a $out/lib
      runHook postInstall
    '';
  };

  gklib = pkgs.stdenv.mkDerivation rec {
    pname = "gklib";
    version = "dev";

    src = pkgs.fetchFromGitHub {
      owner = "KarypisLab";
      repo = "GKlib";
      rev = "6e7951358fd896e2abed7887196b6871aac9f2f8";
      hash = "sha256-jT0hT5Y3E8GnE8OJWzDj5rtz9s59sMEXLduUnBV0I0Y=";
    };

    buildInputs = [
      aux-impls
    ];

    patches = [
      ./normalization/gklib-normalized.patch
    ];

    cmakeFlags = [
      # "-DCMAKE_C_FLAGS_RELEASE=-O2"
      # "-DASSERT=ON" # generally don't want asserts to match release behavior
      "-DNORMALIZED=ON"
    ];

    nativeBuildInputs = with pkgs; [ cmake ];
  };
in
  pkgs.metis.overrideAttrs (finalAttrs: previousAttrs: rec {
    pname = "metis-normalized";
    src = pkgs.fetchFromGitHub {
      owner = "KarypisLab";
      repo = "METIS";
      rev = "v5.2.1";
      hash = "sha256-eddLR6DvZ+2LeR0DkknN6zzRvnW+hLN2qeI+ETUPcac=";
    };

    buildInputs = [
      aux-impls
    ];

    patches = [
      ./normalization/metis-normalized.patch
    ];

    preConfigure = ''
      mkdir -p build/xinclude
      echo "#define IDXTYPEWIDTH 32" > build/xinclude/metis.h
      echo "#define REALTYPEWIDTH 32" >> build/xinclude/metis.h
      cat include/metis.h >> build/xinclude/metis.h
      cp include/CMakeLists.txt build/xinclude
    '';

    cmakeFlags = [
      # "-DCMAKE_C_FLAGS_RELEASE=-O2"
      # "-DCMAKE_C_FLAGS=-w"
      "-DGKLIB_PATH=${gklib}"
      # "-DASSERT=ON" # generally don't want asserts to match release behavior
      "-DNORMALIZED=ON"
      # "-DASSERT2=1"
    ];
  })
