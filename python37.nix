# Nix script to set up appropriate environmen to run python2-based SST sequence analysis code.
#

with import <nixpkgs>{};  # import nixos packages

let 
  nupack=stdenv.mkDerivation{
    name="nupack";
    src=./nupack_viennaRNA/nupack3.0.6.tar.gz;  # nix knows to unzip and untar and run make
    buildPhase=''    # writing a multiline string using 2 pairs of single quotes, that buildPHase intreprets as a bash script
      echo "Running the nupack build phase" 
      make'';
    installPhase=''  # installPhase intreprets double-quoted string as a bash script
      mkdir -p $out; # $out is bash variable that nix sets before running this bash script
      cp -r bin $out;
      cp -r parameters $out;
     '';
  }; 

  cc=stdenv.mkDerivation{
  name="cc";
  phases = [ "installPhase" ];
  installPhase="mkdir -p $out/bin ; cp ${pkgs.gcc.out}/bin/gcc $out/bin/cc";
  };

  viennaRNA=stdenv.mkDerivation{
    name="viennaRNA";
    src=./nupack_viennaRNA/ViennaRNA-2.1.9.tar.gz;
    buildInputs=[which perl gnumake];
    installPhase='' 
    make install
   	'';
  }; in  

stdenv.mkDerivation{
  name="nupack_viennaRNA";
  buildInputs=[coreutils python37 python37Packages.numpy nupack viennaRNA python37Packages.matplotlib]; #       # here I put a list of names of packages I want to install
  NUPACKHOME="${nupack.out}";
}



