
# DNA sequence designer for algorithmic self-assembly of iterated Boolean circuits 

The code in this directory is shipped with the manuscript:
"Diverse and robust molecular algorithms using reprogrammable DNA self-assembly".
Woods*, Doty*, Myhrvold, Hui, Zhou, Yin and Winfree. (*joint first co-authors)
Nature 567:366–372, 2019.

This software can be installed on Linux and MacOS machines in several ways. (Installation on Windows is not possible as a piece of software we require, namely NUPACK, is not supported on Windows.)

We strongly recommend following the Nix instructions outlined below in part A as this provides an automated install procedure that essentially creates an isolated environment with the required prerequisites and without otherwise impacting your system. The nix-installed packages can be easily removed at a later time.

We also provide alternative do-it-yourself instructions (further below in part B) for those who do not wish to install the Nix package manager, but are happy to alter the global state of their system (however, configuration via B might be trickier than A).

Example invocations of the sequence designer are given in part C.

Example usage of our sequence analysis code is given in part D. 

A note on NUPACK and ViennaRNA versions: Our DNA sequence designer used nupack3.0.4 to design the 6bit IBC DNA sequences in the publication, but we include nupack3.0.6 here. We believe that choosing between either version should not affect sequence quality. We have not updated to later versions of NUPACK than nupack3.0.6 because they would require some internal changes to our code to accommodate interface and I/O differences. We used ViennaRNA-2.1.9 in the paper, the same version used here. 

The more recent [nuad](https://github.com/UC-Davis-molecular-computing/nuad/) DNA sequence library makes use of the sequence design principles developed here, and although it includes a uniquely-addressed SST canvas example, as of 2023 it does not include a designer for algorithmic SST self-assembly with the biophysical criteria applied here. 


## A. Nix Instructions

1. Install the [Nix](https://nixos.org/nix/) package manager by following the instructions at <https://nixos.org/download>.
2. During the install, follow any instructions given on screen.
3. Create a directory and place the sequence design code there.
```
git clone https://github.com/DNA-and-Natural-Algorithms-Group/SST-sequence-designer
``` 
4. That directory needs to contain the file default.nix shipped with the sequence design code. You'll need an internet connection. In that directory run the command:
```
nix-shell
```
This will download and install a number of dependencies, essentially in an isolated environment (i.e. not on your system path), including nupack, ViennaRNA, python3 and others. Nix should then present you a new bash shell ready for use. The dependencies required for our sequence designer are locally available to the nix/bash shell and will not be globally available on your system (nor will they interfere with your current system setup).

Running `nix-shell` runs `default.nix`, which uses modern versions of nix. In order to use older versions of nix (2019, or 2023), run: `nix-shell default-2023.nix` and `nix-shell default-2019.nix`. 


## B. Do-it-yourself Instructions

The sequence designer is written in python and relies on the specific point versions of NUPACK and ViennaRNA; using different versions may cause an error. We recommend the sequence designer by run using Python 3 (but also works with Python 2). It requires a number of python packages including numpy and matplotlib.

1. Install python3 (or python2).

2. Install the python packages for numpy and matplotlib.

3. NUPACK: Install NUPACK version 3.0.6. At the time of writing, the source is publicly available by request from the authors here: http://www.nupack.org/downloads (shipped as nupack3.0.6.tar.gz). Extract and run GNU make. The pdf file doc/NUPACK-UserGuide-3.0.pdf included has helpful info on NUPACK. The unix path must be able to find NUPACK executables (e.g. pfunc, mfe), and the environment path variable NUPACKHOME must be set, as described in the NUPACK installation instructions NUPACK-UserGuide-3.0.pdf. The following URL has a mirror of the nupack source:
http://www.dna.caltech.edu/SupplementaryMaterial/Algorithmic_SST/archived_software/nupack3.0.6.tar.gz

4. ViennaRNA: Install ViennaRNA-2.1.9. At the time of writing the source is publicly available here: https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_1_x/ViennaRNA-2.1.9.tar.gz. Follow the instructions in the "INSTALL" text file. You will need to set the environment path variable VIENNARNA_PARAMS_PATH to the location of the ViennaRNA parameters (dna_mathews1999.par and dna_mathews2004.par) which on some systems can be done via (and should be made part of your shell, e.g. bash, startup routine):
export VIENNARNA_PARAMS_PATH=/usr/local/share/ViennaRNA/
The following URL has a mirror of the ViennaRNA sources:
http://www.dna.caltech.edu/SupplementaryMaterial/Algorithmic_SST/archived_software/ViennaRNA-2.1.9.tar.gz


## C. Designing sequences -- example commands

An example command to quickly design a small demo sequence set:
```
python3 atam2ssts.py  -p input_tilesets/demo.py  -o  output_DNA_sequences/demo_sequences.txt
```
Note that the -p (input tile set and energy thresholds) and -o (output sequence file) parameters are mandatory. 

An example command to quickly design 355 sequences using loose parameters that will presumably result in worse sequences than those used in our experimental work:
```
python3 atam2ssts.py  -p input_tilesets/IBC_6bit_loose_params.py  -o output_DNA_sequences/IBC_6bit_loose_params_sequences.txt
```

An example command to design 355 sequences using the parameters used in our experimental work (takes a very long time and may not halt for some runs, depending on random number seed):
```
python3 atam2ssts.py  -p input_tilesets/IBC_6bit.py  -o output_DNA_sequences/IBC_6bit_sequences.txt
```

For help, run:
```
python3 atam2ssts.py --help 
```


## D. Analysing designed sequence -- example commands

Input can be taken directly from a run of sequence designer (e.g. demo_sequences.txt in the example above), or can be an idt-formatted file. The code produces a set of pdf plots that describe thermodynamic properties of the input DNA sequences. 

An example command to analyse a small demo sequence set, at a default temperature of 53.0 C: 
```
python3 analyse_seqs.py output_DNA_sequences/demo_sequences.txt
```
which generates plots that are placed in /output_DNA_sequences/53.0 C/

An example command to analyse a small demo sequence set, at a temperature of 20.0 C: 
```
python3 analyse_seqs.py output_DNA_sequences/demo_sequences.txt -T 20.0
```

For help, run:
```
python3 analyse_seqs.py --help 
```


## E. License

The license that applies to all of the code in this repository is the MIT license see the file "LICENSE". We have also included the NUPACK license (as LICENSE_NUPACK) and ViennaRNA license (as LICENSE_viennaRNA), our python code calls the executable binaries of these software packages but does not link to them directly.



