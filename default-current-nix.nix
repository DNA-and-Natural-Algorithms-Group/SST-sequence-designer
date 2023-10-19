
# Nix script to set up appropriate environment to run python-based SST sequence designer, defaulting to python3. 
# Type 'nix-shell' to run this script. 
# 
# Shipped with DNA single-stranded tile (SST) sequence designer used in the following publication.
# "Diverse and robust molecular algorithms using reprogrammable DNA self-assembly"
# Woods*, Doty*, Myhrvold, Hui, Zhou, Yin, Winfree. (*Joint first co-authors). Nature, 2019

with import <nixpkgs>{};  # import nixos packages

let 
  nupack=stdenv.mkDerivation{
    name="nupack";
    #src=./nupack_viennaRNA/nupack3.0.6.tar.gz;  # nix knows to unzip and untar and run make
    src = fetchurl {
    url = "http://www.dna.caltech.edu/SupplementaryMaterial/Algorithmic_SST/archived_software/nupack3.0.6.tar.gz";
    sha256 = "1kv5irz8n57875dgzr5w3zc9xsy8dbvcfk3iszw5ask942dpdpvw";
    };
    buildPhase = ''
        # -fcommon is required for modern versions of gcc to compile nupack3.0.6 without a 
        # multiple definitions error.
        export NUPACK_CFLAGS="-std=c99 -O3 -Wall -Wmissing-prototypes -Wmissing-declarations -fcommon"
        export NUPACK_CXXFLAGS="-Wall -Wmissing-prototypes -Wmissing-declarations -fcommon"
        CFLAGS=-fcommon make'';
    installPhase = ''  # installPhase interprets double-quoted string as a bash script
          mkdir -p $out; # $out is bash variable that nix sets before running this bash script
          cp -r bin $out;
          cp -r parameters $out;
    '';
  };

  viennaRNA=stdenv.mkDerivation{
    name="viennaRNA";
    src = fetchurl {
    url = "http://www.dna.caltech.edu/SupplementaryMaterial/Algorithmic_SST/archived_software/ViennaRNA-2.1.9.tar.gz";
    sha256 = "1swjnfir5gx424srsnggw4sf8x0p8kiqfzgzp5m34zdzvn4nlzrn";
    };
    # if one has a local copy of the source one can replace the above src instruction with: 
    # src=./ViennaRNA-2.1.9.tar.gz;
    buildInputs = with pkgs; [ which perl gnumake ];
    NIX_CFLAGS_COMPILE = "-fcommon";
    NIX_CXXFLAGS_COMPILE = "-std=c++14";
    buildPhase = ''
      echo "Running the viennaRNA build phase" 
      export CFLAGS="-fcommon"
      export CXXFLAGS="-std=c++14"
      ./configure --prefix=$out --without-perl --without-doc --without-doc-html --without-doc-pdf --without-kinfold --without-svm --without-forester
      make'';
    installPhase = ''
       make install
      	'';
  }; in  

stdenv.mkDerivation{
  name="nupack_viennaRNA";
  buildInputs = with pkgs; [
    coreutils
    (python3.withPackages (ps: [ ps.numpy ps.matplotlib ]))
    nupack
    viennaRNA
  ]; # here I put a list of names of packages I want to install
  NUPACKHOME = "${nupack.out}";
  VIENNARNA_PARAMS_PATH="${viennaRNA.out}/share/ViennaRNA/";
}
