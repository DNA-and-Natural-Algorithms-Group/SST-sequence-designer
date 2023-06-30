# Nix script to set up appropriate environment to run python2-based SST sequence designer.

with import <nixpkgs> { }; # import nixos packages

let
  pkgs = import (builtins.fetchGit {
    name = "2019-03-revision";
    url = "https://github.com/NixOS/nixpkgs/";
    ref = "refs/heads/nixpkgs-unstable";
    rev = "25b53f32236b308172673c958f570b5b488c7b73";
  }) { };
  nupack = stdenv.mkDerivation {
    name = "nupack";
    src =
      ./nupack_viennaRNA/nupack3.0.6.tar.gz; # nix knows to unzip and untar and run make
    #postPatch = ''
    #substituteInPlace src/thermo/utils/init.c --replace "getenv(\"NUPACKHOME\")" "\"$out\""
    #substituteInPlace src/thermo/distributions/ReadCommandLine.c --replace "getenv(\"NUPACKHOME\")" "\"$out\""
    #'';
    buildPhase = ''
      # writing a multiline string using 2 pairs of single quotes, that buildPHase intreprets as a bash script
        echo "Running the nupack build phase" 

        # -fcommon is required for modern versions of gcc to compile nupack3.0.6 without a 
        # multiple definitions error.
        export NUPACK_CFLAGS="-std=c99 -O3 -Wall -Wmissing-prototypes -Wmissing-declarations -fcommon"
        export NUPACK_CXXFLAGS="-Wall -Wmissing-prototypes -Wmissing-declarations -fcommon"
        CFLAGS=-fcommon make'';
    installPhase = ''
      # installPhase intreprets double-quoted string as a bash script
          mkdir -p $out; # $out is bash variable that nix sets before running this bash script
          cp -r bin $out;
          cp -r parameters $out;
    '';
  };

  viennaRNA = stdenv.mkDerivation {
    name = "viennaRNA";
    src = ./nupack_viennaRNA/ViennaRNA-2.1.9.tar.gz;
    buildInputs = [ which perl gnumake ];
    NIX_CFLAGS_COMPILE = "-fcommon";
    NIX_CXXFLAGS_COMPILE = "-std=c++14";
    buildPhase = ''
      echo "Running the viennaRNA build phase" 
      export CFLAGS="-fcommon"
      export CXXFLAGS="-std=c++14"
      ./configure --prefix=$out
      make'';
    installPhase = ''
       make install
      	'';
  };

in stdenv.mkDerivation {
  name = "nupack_viennaRNA";
  buildInputs = with pkgs; [
    coreutils
    (python3.withPackages (ps: [ ps.numpy ps.matplotlib ]))
    nupack
    viennaRNA
  ]; # here I put a list of names of packages I want to install
  NUPACKHOME = "${nupack.out}";
}

# for analysis  python2Packages.matplotlib

