# Nix script to set up appropriate environmen to run python2-based SST sequence designer.
#
# 
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
     '';
  }; 

  cc=stdenv.mkDerivation{
  name="cc";
  phases = [ "installPhase" ];
  installPhase="mkdir -p $out/bin ; cp ${pkgs.gcc.out}/bin/gcc $out/bin/cc";
  };

  viennaRNA=stdenv.mkDerivation{
    name="viennaRNA";
    #preConfigure="export CC=${pkgs.gcc.out}/bin/gcc";
    #preBuild="export CC=${pkgs.gcc.out}/bin/gcc";
    src=./nupack_viennaRNA/ViennaRNA-2.1.9.tar.gz;
    buildInputs=[which perl gnumake];
    buildPhase=''    # writing a multiline string using 2 pairs of single quotes, that buildPHase intreprets as a bash script
      echo "Running the viennaRNA build phase" 
      #./configure
      #make
      #make check
      '';

    installPhase='' 
    make install
    #make clean
    # The following bash commands copy the compiled binaries to the output.
    # We found (by running "find ." here) that the binaries are compiled in the 
    # build stage and placed in the directory ./src/bin . Different systems provide
    # different versions of ls, so we include coreutils below to get whatever is 
    # the ls shipped with nix coreutils 
    #mkdir -p $out/bin 
    #for i in $(ls src/bin/RNA* | grep -v "\." ); do
    #  cp $i $out/bin
    #done 
   	'';
  }; in  

stdenv.mkDerivation{
  name="nupack_viennaRNA";
  buildInputs=[coreutils python2 python2Packages.numpy nupack viennaRNA];   # here I put a list of names of packages I want to install
}



