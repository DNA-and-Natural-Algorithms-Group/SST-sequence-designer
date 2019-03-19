
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

    installPhase=''  # installPhase interprets double-quoted string as a bash script
      mkdir -p $out; # $out is a bash variable that nix sets before running this bash script
      cp -r bin $out;
      cp -r parameters $out;
     '';
  }; 

  #cc=stdenv.mkDerivation{
  #name="cc";
  #phases = [ "installPhase" ];
  #installPhase="mkdir -p $out/bin ; cp ${pkgs.gcc.out}/bin/gcc $out/bin/cc";
  #};

  viennaRNA=stdenv.mkDerivation{
    name="viennaRNA";
    src = fetchurl {
    url = "http://www.dna.caltech.edu/SupplementaryMaterial/Algorithmic_SST/archived_software/ViennaRNA-2.1.9.tar.gz";
    sha256 = "1swjnfir5gx424srsnggw4sf8x0p8kiqfzgzp5m34zdzvn4nlzrn";
    };
    # if one has a local copy of the source one can replace the above src instruction with: 
    # src=./ViennaRNA-2.1.9.tar.gz;

  buildInputs=[which perl gnumake];

  }; in  


stdenv.mkDerivation{
  name="nupack_viennaRNA";
  # python3:
  buildInputs=[coreutils python3 python3Packages.numpy python3Packages.matplotlib nupack viennaRNA];
  # python2:
  # buildInputs=[coreutils python2 python2Packages.numpy python2Packages.matplotlib nupack viennaRNA]; 
  NUPACKHOME="${nupack.out}";
  VIENNARNA_PARAMS_PATH="${viennaRNA.out}/share/ViennaRNA/";
}

