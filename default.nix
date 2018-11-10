# Nix script to set up appropriate environmen to run python2-based SST sequence designer.
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
  #buildInputs=[coreutils python2 python2Packages.numpy nupack viennaRNA]; #   python2Packages.matplotlib python3Packages.matplotlib  python3 python3Packages.numpy   # here I put a list of names of packages I want to install
  buildInputs=[coreutils python3 python3Packages.numpy nupack viennaRNA python3Packages.matplotlib]; 
  NUPACKHOME="${nupack.out}";
}



