
{ pkgs ? import <nixpkgs> {} }:
pkgs.stdenv.mkDerivation rec {
  pname = "aux-normalization";
  version = "dev";

  src = ./.;

  nativeBuildInputs = with pkgs; [ rustc ];

  buildPhase = ''
      runHook preBuild
      rustc functions.rs
      runHook postBuild
      '';

  installPhase = ''
      runHook preInstall
      mkdir -p $out/lib
      mv libfunctions.so $out/lib
      ln -s $out/lib/libfunctions.so $out/lib/libnormalization-functions.so

      runHook postInstall
      '';
}
