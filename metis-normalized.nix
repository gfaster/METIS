{ 
aux-normalization,
pkgs ? import <nixpkgs> {}
}:
let
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
      aux-normalization
    ];

    patches = [
      ./normalization/gklib-normalized.patch
    ];

    cmakeFlags = [
      "-DCMAKE_C_FLAGS_RELEASE=-O1"
      # "-DASSERT=ON" # generally don't want asserts to match release behavior
      "-DNORMALIZED=ON"
    ];

    nativeBuildInputs = with pkgs; [ cmake ];
  };

  metis = pkgs.stdenv.mkDerivation rec {
    pname = "metis-normalized";
    version = "5.2.1";

    src = pkgs.fetchFromGitHub {
      owner = "KarypisLab";
      repo = "METIS";
      rev = "a6e6a2cfa92f93a3ee2971ebc9ddfc3b0b581ab2";
      hash = "sha256-VObeibIQtZ8bFnFo2jXsCh8pqEdVkVOtCBBSvZ4lDEM=";
    };

    nativeBuldInputs = [ ];

    # for whatever reason, cmake has to be in buildInputs and doesn't work in
    # nativeBuildInputs
    buildInputs = [
      pkgs.cmake
      aux-normalization
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

    dontStrip = true;

    cmakeFlags = [
      "-DCMAKE_C_FLAGS_RELEASE=-O1"
      "-DCMAKE_C_FLAGS=-w"
      "-DGKLIB_PATH=${gklib}"
      # "-DASSERT=ON" # generally don't want asserts to match release behavior
      "-DNORMALIZED=ON"
      # "-DASSERT2=ON"
    ];
    };
in
  metis
  # pkgs.enableDebugging (pkgs.callPackage metis {})
