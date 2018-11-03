
# Input file that includes input tileset and threshold parameters for SST sequence designer. 
#
# Specification of 6-bit IBC tile set with 2x2 proofreading with energy thresholds 
# as used in the following publication.
#
# Woods*, Doty*, Myhrvold, Hui, Zhou, Yin, Winfree. (*Joint first co-authors)
# Diverse and robust molecular algorithms using reprogrammable DNA self-assembly
# 
# Note that all energy parameter values are sign-flipped (i.e. more positive is more favourable)



# "global" energy interval for all glues, except those glues with user-specified sequences;
# evaluated using a simple nearest-neighbour model
lowDG=8.9
highDG=9.2
biotin_strength_boost=1.1

# temperature in degrees C at which to evaluate energies
temperature = 53.0
    
# allowed energy of secondary structure of each individual tile (via nupack)
tile_sec_struct = 1.65 

# allowed energy of secondary structure of each pair of tiles that don't share any glues (via RNAduplex)
tile_orth = 5.4  

# allowed energy of secondary structure of pair of tiles that share glues (via RNAduplex)
tile_orth_share = 7.4 

# threshold for allowed energy of interaction between algorithmically conflicting glues
# that might be co-located during a strength-1 binding event on an otherwise perfect lattice (via nupack)
orth_algorithmic_conflict = 1.6

# threshold for allowed energy of interaction between algorithmically conflicting glues
# that might be co-located during a strength-1 binding event on a lattice that possibly 
# has some other kinds of error (via nupack)
orth_algorithmic_conflict_generalized = 2.6

# threshold for allowed energy of glues that will be colocated during correct growth 
# (must be broken up to allow binding of new tile) (via nupack)
orth_colocated = 1.4

# tiles should bind to lattice with energy strictly greater than the following thresold (via nupack)
lattice_binding_lower_threshold = 12.3


# if this line is present, then parity 0 tiles will use a 3-letter code
# A,C,T, and parity 1 tiles will use A,G,T
three_letter_code = True

# number of bases that are exceptions to the three letter code, e.g. TGACTGGTAAAG has 
# three letter code A,T,G with a three_letter_code_exceptions = 1
three_letter_code_exceptions = 1

# if endGC == True then each domain begins and end with G or C, e.g. GTTAAGCGTCC 
endGC = False
# if endAT == True then each domain begins and end with A or T, e.g. TTTAAGCGTCA
endAT = True

# when designing deterministic tilesets setting detect_nondeterminism=True performs 
# a syntax check on the tile set to ensure there is no nondeterminism 
detect_nondeterminism=False



# List of tiles, each tile is a python dictionary
tiles = \
[{'name':'U6;00->00;nw' , 'parity':1, 'S':'iL:U6;00->00' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;00->00' },
 {'name':'U6;00->00;ne' , 'parity':0, 'S':'iR:U6;00->00' , 'W':'iT:U6;00->00' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_0'     },
 {'name':'U6;00->00;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;00->00' , 'E':'iB:U6;00->00' },
 {'name':'U6;00->00;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;00->00' , 'N':'iR:U6;00->00' , 'E':'B:U6_E_0'     },
 {'name':'U6;00->01;nw' , 'parity':1, 'S':'iL:U6;00->01' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;00->01' },
 {'name':'U6;00->01;ne' , 'parity':0, 'S':'iR:U6;00->01' , 'W':'iT:U6;00->01' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_1'     },
 {'name':'U6;00->01;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;00->01' , 'E':'iB:U6;00->01' },
 {'name':'U6;00->01;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;00->01' , 'N':'iR:U6;00->01' , 'E':'B:U6_E_1'     },
 {'name':'U6;00->10;nw' , 'parity':1, 'S':'iL:U6;00->10' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;00->10' },
 {'name':'U6;00->10;ne' , 'parity':0, 'S':'iR:U6;00->10' , 'W':'iT:U6;00->10' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_0'     },
 {'name':'U6;00->10;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;00->10' , 'E':'iB:U6;00->10' },
 {'name':'U6;00->10;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;00->10' , 'N':'iR:U6;00->10' , 'E':'B:U6_E_0'     },
 {'name':'U6;00->11;nw' , 'parity':1, 'S':'iL:U6;00->11' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;00->11' },
 {'name':'U6;00->11;ne' , 'parity':0, 'S':'iR:U6;00->11' , 'W':'iT:U6;00->11' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_1'     },
 {'name':'U6;00->11;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;00->11' , 'E':'iB:U6;00->11' },
 {'name':'U6;00->11;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;00->11' , 'N':'iR:U6;00->11' , 'E':'B:U6_E_1'     },
 {'name':'U6;01->00;nw' , 'parity':1, 'S':'iL:U6;01->00' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;01->00' },
 {'name':'U6;01->00;ne' , 'parity':0, 'S':'iR:U6;01->00' , 'W':'iT:U6;01->00' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_0'     },
 {'name':'U6;01->00;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;01->00' , 'E':'iB:U6;01->00' },
 {'name':'U6;01->00;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;01->00' , 'N':'iR:U6;01->00' , 'E':'B:U6_E_0'     },
 {'name':'U6;01->01;nw' , 'parity':1, 'S':'iL:U6;01->01' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;01->01' },
 {'name':'U6;01->01;ne' , 'parity':0, 'S':'iR:U6;01->01' , 'W':'iT:U6;01->01' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_1'     },
 {'name':'U6;01->01;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;01->01' , 'E':'iB:U6;01->01' },
 {'name':'U6;01->01;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;01->01' , 'N':'iR:U6;01->01' , 'E':'B:U6_E_1'     },
 {'name':'U6;01->10;nw' , 'parity':1, 'S':'iL:U6;01->10' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;01->10' },
 {'name':'U6;01->10;ne' , 'parity':0, 'S':'iR:U6;01->10' , 'W':'iT:U6;01->10' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_0'     },
 {'name':'U6;01->10;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;01->10' , 'E':'iB:U6;01->10' },
 {'name':'U6;01->10;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;01->10' , 'N':'iR:U6;01->10' , 'E':'B:U6_E_0'     },
 {'name':'U6;01->11;nw' , 'parity':1, 'S':'iL:U6;01->11' , 'W':'T:U6_W_0'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;01->11' },
 {'name':'U6;01->11;ne' , 'parity':0, 'S':'iR:U6;01->11' , 'W':'iT:U6;01->11' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_1'     },
 {'name':'U6;01->11;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_0'     , 'N':'iL:U6;01->11' , 'E':'iB:U6;01->11' },
 {'name':'U6;01->11;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;01->11' , 'N':'iR:U6;01->11' , 'E':'B:U6_E_1'     },
 {'name':'U6;10->00;nw' , 'parity':1, 'S':'iL:U6;10->00' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;10->00' },
 {'name':'U6;10->00;ne' , 'parity':0, 'S':'iR:U6;10->00' , 'W':'iT:U6;10->00' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_0'     },
 {'name':'U6;10->00;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;10->00' , 'E':'iB:U6;10->00' },
 {'name':'U6;10->00;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;10->00' , 'N':'iR:U6;10->00' , 'E':'B:U6_E_0'     },
 {'name':'U6;10->01;nw' , 'parity':1, 'S':'iL:U6;10->01' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;10->01' },
 {'name':'U6;10->01;ne' , 'parity':0, 'S':'iR:U6;10->01' , 'W':'iT:U6;10->01' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_1'     },
 {'name':'U6;10->01;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;10->01' , 'E':'iB:U6;10->01' },
 {'name':'U6;10->01;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;10->01' , 'N':'iR:U6;10->01' , 'E':'B:U6_E_1'     },
 {'name':'U6;10->10;nw' , 'parity':1, 'S':'iL:U6;10->10' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;10->10' },
 {'name':'U6;10->10;ne' , 'parity':0, 'S':'iR:U6;10->10' , 'W':'iT:U6;10->10' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_0'     },
 {'name':'U6;10->10;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;10->10' , 'E':'iB:U6;10->10' },
 {'name':'U6;10->10;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;10->10' , 'N':'iR:U6;10->10' , 'E':'B:U6_E_0'     },
 {'name':'U6;10->11;nw' , 'parity':1, 'S':'iL:U6;10->11' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_0'     , 'E':'iT:U6;10->11' },
 {'name':'U6;10->11;ne' , 'parity':0, 'S':'iR:U6;10->11' , 'W':'iT:U6;10->11' , 'N':'R:U6_N_0'     , 'E':'T:U6_E_1'     },
 {'name':'U6;10->11;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;10->11' , 'E':'iB:U6;10->11' },
 {'name':'U6;10->11;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;10->11' , 'N':'iR:U6;10->11' , 'E':'B:U6_E_1'     },
 {'name':'U6;11->00;nw' , 'parity':1, 'S':'iL:U6;11->00' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;11->00' },
 {'name':'U6;11->00;ne' , 'parity':0, 'S':'iR:U6;11->00' , 'W':'iT:U6;11->00' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_0'     },
 {'name':'U6;11->00;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;11->00' , 'E':'iB:U6;11->00' },
 {'name':'U6;11->00;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;11->00' , 'N':'iR:U6;11->00' , 'E':'B:U6_E_0'     },
 {'name':'U6;11->01;nw' , 'parity':1, 'S':'iL:U6;11->01' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;11->01' },
 {'name':'U6;11->01;ne' , 'parity':0, 'S':'iR:U6;11->01' , 'W':'iT:U6;11->01' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_1'     },
 {'name':'U6;11->01;sw' , 'parity':0, 'S':'L:U6_S_0'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;11->01' , 'E':'iB:U6;11->01' },
 {'name':'U6;11->01;se' , 'parity':1, 'S':'R:U6_S_0'     , 'W':'iB:U6;11->01' , 'N':'iR:U6;11->01' , 'E':'B:U6_E_1'     },
 {'name':'U6;11->10;nw' , 'parity':1, 'S':'iL:U6;11->10' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;11->10' },
 {'name':'U6;11->10;ne' , 'parity':0, 'S':'iR:U6;11->10' , 'W':'iT:U6;11->10' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_0'     },
 {'name':'U6;11->10;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;11->10' , 'E':'iB:U6;11->10' },
 {'name':'U6;11->10;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;11->10' , 'N':'iR:U6;11->10' , 'E':'B:U6_E_0'     },
 {'name':'U6;11->11;nw' , 'parity':1, 'S':'iL:U6;11->11' , 'W':'T:U6_W_1'     , 'N':'L:U6_N_1'     , 'E':'iT:U6;11->11' },
 {'name':'U6;11->11;ne' , 'parity':0, 'S':'iR:U6;11->11' , 'W':'iT:U6;11->11' , 'N':'R:U6_N_1'     , 'E':'T:U6_E_1'     },
 {'name':'U6;11->11;sw' , 'parity':0, 'S':'L:U6_S_1'     , 'W':'B:U6_W_1'     , 'N':'iL:U6;11->11' , 'E':'iB:U6;11->11' },
 {'name':'U6;11->11;se' , 'parity':1, 'S':'R:U6_S_1'     , 'W':'iB:U6;11->11' , 'N':'iR:U6;11->11' , 'E':'B:U6_E_1'     },
 {'name':'U4;00->00;nw' , 'parity':1, 'S':'iL:U4;00->00' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;00->00' },
 {'name':'U4;00->00;ne' , 'parity':0, 'S':'iR:U4;00->00' , 'W':'iT:U4;00->00' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_0'     },
 {'name':'U4;00->00;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;00->00' , 'E':'iB:U4;00->00' },
 {'name':'U4;00->00;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;00->00' , 'N':'iR:U4;00->00' , 'E':'B:U4_E_0'     },
 {'name':'U4;00->01;nw' , 'parity':1, 'S':'iL:U4;00->01' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;00->01' },
 {'name':'U4;00->01;ne' , 'parity':0, 'S':'iR:U4;00->01' , 'W':'iT:U4;00->01' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_1'     },
 {'name':'U4;00->01;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;00->01' , 'E':'iB:U4;00->01' },
 {'name':'U4;00->01;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;00->01' , 'N':'iR:U4;00->01' , 'E':'B:U4_E_1'     },
 {'name':'U4;00->10;nw' , 'parity':1, 'S':'iL:U4;00->10' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;00->10' },
 {'name':'U4;00->10;ne' , 'parity':0, 'S':'iR:U4;00->10' , 'W':'iT:U4;00->10' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_0'     },
 {'name':'U4;00->10;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;00->10' , 'E':'iB:U4;00->10' },
 {'name':'U4;00->10;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;00->10' , 'N':'iR:U4;00->10' , 'E':'B:U4_E_0'     },
 {'name':'U4;00->11;nw' , 'parity':1, 'S':'iL:U4;00->11' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;00->11' },
 {'name':'U4;00->11;ne' , 'parity':0, 'S':'iR:U4;00->11' , 'W':'iT:U4;00->11' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_1'     },
 {'name':'U4;00->11;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;00->11' , 'E':'iB:U4;00->11' },
 {'name':'U4;00->11;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;00->11' , 'N':'iR:U4;00->11' , 'E':'B:U4_E_1'     },
 {'name':'U4;01->00;nw' , 'parity':1, 'S':'iL:U4;01->00' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;01->00' },
 {'name':'U4;01->00;ne' , 'parity':0, 'S':'iR:U4;01->00' , 'W':'iT:U4;01->00' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_0'     },
 {'name':'U4;01->00;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;01->00' , 'E':'iB:U4;01->00' },
 {'name':'U4;01->00;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;01->00' , 'N':'iR:U4;01->00' , 'E':'B:U4_E_0'     },
 {'name':'U4;01->01;nw' , 'parity':1, 'S':'iL:U4;01->01' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;01->01' },
 {'name':'U4;01->01;ne' , 'parity':0, 'S':'iR:U4;01->01' , 'W':'iT:U4;01->01' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_1'     },
 {'name':'U4;01->01;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;01->01' , 'E':'iB:U4;01->01' },
 {'name':'U4;01->01;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;01->01' , 'N':'iR:U4;01->01' , 'E':'B:U4_E_1'     },
 {'name':'U4;01->10;nw' , 'parity':1, 'S':'iL:U4;01->10' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;01->10' },
 {'name':'U4;01->10;ne' , 'parity':0, 'S':'iR:U4;01->10' , 'W':'iT:U4;01->10' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_0'     },
 {'name':'U4;01->10;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;01->10' , 'E':'iB:U4;01->10' },
 {'name':'U4;01->10;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;01->10' , 'N':'iR:U4;01->10' , 'E':'B:U4_E_0'     },
 {'name':'U4;01->11;nw' , 'parity':1, 'S':'iL:U4;01->11' , 'W':'T:U4_W_0'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;01->11' },
 {'name':'U4;01->11;ne' , 'parity':0, 'S':'iR:U4;01->11' , 'W':'iT:U4;01->11' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_1'     },
 {'name':'U4;01->11;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_0'     , 'N':'iL:U4;01->11' , 'E':'iB:U4;01->11' },
 {'name':'U4;01->11;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;01->11' , 'N':'iR:U4;01->11' , 'E':'B:U4_E_1'     },
 {'name':'U4;10->00;nw' , 'parity':1, 'S':'iL:U4;10->00' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;10->00' },
 {'name':'U4;10->00;ne' , 'parity':0, 'S':'iR:U4;10->00' , 'W':'iT:U4;10->00' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_0'     },
 {'name':'U4;10->00;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;10->00' , 'E':'iB:U4;10->00' },
 {'name':'U4;10->00;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;10->00' , 'N':'iR:U4;10->00' , 'E':'B:U4_E_0'     },
 {'name':'U4;10->01;nw' , 'parity':1, 'S':'iL:U4;10->01' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;10->01' },
 {'name':'U4;10->01;ne' , 'parity':0, 'S':'iR:U4;10->01' , 'W':'iT:U4;10->01' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_1'     },
 {'name':'U4;10->01;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;10->01' , 'E':'iB:U4;10->01' },
 {'name':'U4;10->01;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;10->01' , 'N':'iR:U4;10->01' , 'E':'B:U4_E_1'     },
 {'name':'U4;10->10;nw' , 'parity':1, 'S':'iL:U4;10->10' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;10->10' },
 {'name':'U4;10->10;ne' , 'parity':0, 'S':'iR:U4;10->10' , 'W':'iT:U4;10->10' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_0'     },
 {'name':'U4;10->10;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;10->10' , 'E':'iB:U4;10->10' },
 {'name':'U4;10->10;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;10->10' , 'N':'iR:U4;10->10' , 'E':'B:U4_E_0'     },
 {'name':'U4;10->11;nw' , 'parity':1, 'S':'iL:U4;10->11' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_0'     , 'E':'iT:U4;10->11' },
 {'name':'U4;10->11;ne' , 'parity':0, 'S':'iR:U4;10->11' , 'W':'iT:U4;10->11' , 'N':'R:U4_N_0'     , 'E':'T:U4_E_1'     },
 {'name':'U4;10->11;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;10->11' , 'E':'iB:U4;10->11' },
 {'name':'U4;10->11;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;10->11' , 'N':'iR:U4;10->11' , 'E':'B:U4_E_1'     },
 {'name':'U4;11->00;nw' , 'parity':1, 'S':'iL:U4;11->00' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;11->00' },
 {'name':'U4;11->00;ne' , 'parity':0, 'S':'iR:U4;11->00' , 'W':'iT:U4;11->00' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_0'     },
 {'name':'U4;11->00;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;11->00' , 'E':'iB:U4;11->00' },
 {'name':'U4;11->00;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;11->00' , 'N':'iR:U4;11->00' , 'E':'B:U4_E_0'     },
 {'name':'U4;11->01;nw' , 'parity':1, 'S':'iL:U4;11->01' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;11->01' },
 {'name':'U4;11->01;ne' , 'parity':0, 'S':'iR:U4;11->01' , 'W':'iT:U4;11->01' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_1'     },
 {'name':'U4;11->01;sw' , 'parity':0, 'S':'L:U4_S_0'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;11->01' , 'E':'iB:U4;11->01' },
 {'name':'U4;11->01;se' , 'parity':1, 'S':'R:U4_S_0'     , 'W':'iB:U4;11->01' , 'N':'iR:U4;11->01' , 'E':'B:U4_E_1'     },
 {'name':'U4;11->10;nw' , 'parity':1, 'S':'iL:U4;11->10' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;11->10' },
 {'name':'U4;11->10;ne' , 'parity':0, 'S':'iR:U4;11->10' , 'W':'iT:U4;11->10' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_0'     },
 {'name':'U4;11->10;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;11->10' , 'E':'iB:U4;11->10' },
 {'name':'U4;11->10;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;11->10' , 'N':'iR:U4;11->10' , 'E':'B:U4_E_0'     },
 {'name':'U4;11->11;nw' , 'parity':1, 'S':'iL:U4;11->11' , 'W':'T:U4_W_1'     , 'N':'L:U4_N_1'     , 'E':'iT:U4;11->11' },
 {'name':'U4;11->11;ne' , 'parity':0, 'S':'iR:U4;11->11' , 'W':'iT:U4;11->11' , 'N':'R:U4_N_1'     , 'E':'T:U4_E_1'     },
 {'name':'U4;11->11;sw' , 'parity':0, 'S':'L:U4_S_1'     , 'W':'B:U4_W_1'     , 'N':'iL:U4;11->11' , 'E':'iB:U4;11->11' },
 {'name':'U4;11->11;se' , 'parity':1, 'S':'R:U4_S_1'     , 'W':'iB:U4;11->11' , 'N':'iR:U4;11->11' , 'E':'B:U4_E_1'     },
 {'name':'U7;00->00;nw' , 'parity':1, 'S':'iL:U7;00->00' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;00->00' },
 {'name':'U7;00->00;ne' , 'parity':0, 'S':'iR:U7;00->00' , 'W':'iT:U7;00->00' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_0'     },
 {'name':'U7;00->00;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;00->00' , 'E':'iB:U7;00->00' },
 {'name':'U7;00->00;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;00->00' , 'N':'iR:U7;00->00' , 'E':'B:U8_W_0'     },
 {'name':'U7;00->01;nw' , 'parity':1, 'S':'iL:U7;00->01' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;00->01' },
 {'name':'U7;00->01;ne' , 'parity':0, 'S':'iR:U7;00->01' , 'W':'iT:U7;00->01' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_1'     },
 {'name':'U7;00->01;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;00->01' , 'E':'iB:U7;00->01' },
 {'name':'U7;00->01;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;00->01' , 'N':'iR:U7;00->01' , 'E':'B:U8_W_1'     },
 {'name':'U7;00->10;nw' , 'parity':1, 'S':'iL:U7;00->10' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;00->10' },
 {'name':'U7;00->10;ne' , 'parity':0, 'S':'iR:U7;00->10' , 'W':'iT:U7;00->10' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_0'     },
 {'name':'U7;00->10;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;00->10' , 'E':'iB:U7;00->10' },
 {'name':'U7;00->10;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;00->10' , 'N':'iR:U7;00->10' , 'E':'B:U8_W_0'     },
 {'name':'U7;00->11;nw' , 'parity':1, 'S':'iL:U7;00->11' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;00->11' },
 {'name':'U7;00->11;ne' , 'parity':0, 'S':'iR:U7;00->11' , 'W':'iT:U7;00->11' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_1'     },
 {'name':'U7;00->11;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;00->11' , 'E':'iB:U7;00->11' },
 {'name':'U7;00->11;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;00->11' , 'N':'iR:U7;00->11' , 'E':'B:U8_W_1'     },
 {'name':'U7;01->00;nw' , 'parity':1, 'S':'iL:U7;01->00' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;01->00' },
 {'name':'U7;01->00;ne' , 'parity':0, 'S':'iR:U7;01->00' , 'W':'iT:U7;01->00' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_0'     },
 {'name':'U7;01->00;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;01->00' , 'E':'iB:U7;01->00' },
 {'name':'U7;01->00;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;01->00' , 'N':'iR:U7;01->00' , 'E':'B:U8_W_0'     },
 {'name':'U7;01->01;nw' , 'parity':1, 'S':'iL:U7;01->01' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;01->01' },
 {'name':'U7;01->01;ne' , 'parity':0, 'S':'iR:U7;01->01' , 'W':'iT:U7;01->01' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_1'     },
 {'name':'U7;01->01;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;01->01' , 'E':'iB:U7;01->01' },
 {'name':'U7;01->01;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;01->01' , 'N':'iR:U7;01->01' , 'E':'B:U8_W_1'     },
 {'name':'U7;01->10;nw' , 'parity':1, 'S':'iL:U7;01->10' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;01->10' },
 {'name':'U7;01->10;ne' , 'parity':0, 'S':'iR:U7;01->10' , 'W':'iT:U7;01->10' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_0'     },
 {'name':'U7;01->10;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;01->10' , 'E':'iB:U7;01->10' },
 {'name':'U7;01->10;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;01->10' , 'N':'iR:U7;01->10' , 'E':'B:U8_W_0'     },
 {'name':'U7;01->11;nw' , 'parity':1, 'S':'iL:U7;01->11' , 'W':'T:U6_E_0'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;01->11' },
 {'name':'U7;01->11;ne' , 'parity':0, 'S':'iR:U7;01->11' , 'W':'iT:U7;01->11' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_1'     },
 {'name':'U7;01->11;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_0'     , 'N':'iL:U7;01->11' , 'E':'iB:U7;01->11' },
 {'name':'U7;01->11;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;01->11' , 'N':'iR:U7;01->11' , 'E':'B:U8_W_1'     },
 {'name':'U7;10->00;nw' , 'parity':1, 'S':'iL:U7;10->00' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;10->00' },
 {'name':'U7;10->00;ne' , 'parity':0, 'S':'iR:U7;10->00' , 'W':'iT:U7;10->00' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_0'     },
 {'name':'U7;10->00;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;10->00' , 'E':'iB:U7;10->00' },
 {'name':'U7;10->00;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;10->00' , 'N':'iR:U7;10->00' , 'E':'B:U8_W_0'     },
 {'name':'U7;10->01;nw' , 'parity':1, 'S':'iL:U7;10->01' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;10->01' },
 {'name':'U7;10->01;ne' , 'parity':0, 'S':'iR:U7;10->01' , 'W':'iT:U7;10->01' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_1'     },
 {'name':'U7;10->01;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;10->01' , 'E':'iB:U7;10->01' },
 {'name':'U7;10->01;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;10->01' , 'N':'iR:U7;10->01' , 'E':'B:U8_W_1'     },
 {'name':'U7;10->10;nw' , 'parity':1, 'S':'iL:U7;10->10' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;10->10' },
 {'name':'U7;10->10;ne' , 'parity':0, 'S':'iR:U7;10->10' , 'W':'iT:U7;10->10' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_0'     },
 {'name':'U7;10->10;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;10->10' , 'E':'iB:U7;10->10' },
 {'name':'U7;10->10;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;10->10' , 'N':'iR:U7;10->10' , 'E':'B:U8_W_0'     },
 {'name':'U7;10->11;nw' , 'parity':1, 'S':'iL:U7;10->11' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_0'     , 'E':'iT:U7;10->11' },
 {'name':'U7;10->11;ne' , 'parity':0, 'S':'iR:U7;10->11' , 'W':'iT:U7;10->11' , 'N':'R:U8_S_0'     , 'E':'T:U8_W_1'     },
 {'name':'U7;10->11;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;10->11' , 'E':'iB:U7;10->11' },
 {'name':'U7;10->11;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;10->11' , 'N':'iR:U7;10->11' , 'E':'B:U8_W_1'     },
 {'name':'U7;11->00;nw' , 'parity':1, 'S':'iL:U7;11->00' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;11->00' },
 {'name':'U7;11->00;ne' , 'parity':0, 'S':'iR:U7;11->00' , 'W':'iT:U7;11->00' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_0'     },
 {'name':'U7;11->00;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;11->00' , 'E':'iB:U7;11->00' },
 {'name':'U7;11->00;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;11->00' , 'N':'iR:U7;11->00' , 'E':'B:U8_W_0'     },
 {'name':'U7;11->01;nw' , 'parity':1, 'S':'iL:U7;11->01' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;11->01' },
 {'name':'U7;11->01;ne' , 'parity':0, 'S':'iR:U7;11->01' , 'W':'iT:U7;11->01' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_1'     },
 {'name':'U7;11->01;sw' , 'parity':0, 'S':'L:U6_N_0'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;11->01' , 'E':'iB:U7;11->01' },
 {'name':'U7;11->01;se' , 'parity':1, 'S':'R:U6_N_0'     , 'W':'iB:U7;11->01' , 'N':'iR:U7;11->01' , 'E':'B:U8_W_1'     },
 {'name':'U7;11->10;nw' , 'parity':1, 'S':'iL:U7;11->10' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;11->10' },
 {'name':'U7;11->10;ne' , 'parity':0, 'S':'iR:U7;11->10' , 'W':'iT:U7;11->10' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_0'     },
 {'name':'U7;11->10;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;11->10' , 'E':'iB:U7;11->10' },
 {'name':'U7;11->10;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;11->10' , 'N':'iR:U7;11->10' , 'E':'B:U8_W_0'     },
 {'name':'U7;11->11;nw' , 'parity':1, 'S':'iL:U7;11->11' , 'W':'T:U6_E_1'     , 'N':'L:U8_S_1'     , 'E':'iT:U7;11->11' },
 {'name':'U7;11->11;ne' , 'parity':0, 'S':'iR:U7;11->11' , 'W':'iT:U7;11->11' , 'N':'R:U8_S_1'     , 'E':'T:U8_W_1'     },
 {'name':'U7;11->11;sw' , 'parity':0, 'S':'L:U6_N_1'     , 'W':'B:U6_E_1'     , 'N':'iL:U7;11->11' , 'E':'iB:U7;11->11' },
 {'name':'U7;11->11;se' , 'parity':1, 'S':'R:U6_N_1'     , 'W':'iB:U7;11->11' , 'N':'iR:U7;11->11' , 'E':'B:U8_W_1'     },
 {'name':'U5;00->00;nw' , 'parity':1, 'S':'iL:U5;00->00' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;00->00' },
 {'name':'U5;00->00;ne' , 'parity':0, 'S':'iR:U5;00->00' , 'W':'iT:U5;00->00' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_0'     },
 {'name':'U5;00->00;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;00->00' , 'E':'iB:U5;00->00' },
 {'name':'U5;00->00;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;00->00' , 'N':'iR:U5;00->00' , 'E':'B:U6_W_0'     },
 {'name':'U5;00->01;nw' , 'parity':1, 'S':'iL:U5;00->01' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;00->01' },
 {'name':'U5;00->01;ne' , 'parity':0, 'S':'iR:U5;00->01' , 'W':'iT:U5;00->01' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_1'     },
 {'name':'U5;00->01;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;00->01' , 'E':'iB:U5;00->01' },
 {'name':'U5;00->01;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;00->01' , 'N':'iR:U5;00->01' , 'E':'B:U6_W_1'     },
 {'name':'U5;00->10;nw' , 'parity':1, 'S':'iL:U5;00->10' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;00->10' },
 {'name':'U5;00->10;ne' , 'parity':0, 'S':'iR:U5;00->10' , 'W':'iT:U5;00->10' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_0'     },
 {'name':'U5;00->10;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;00->10' , 'E':'iB:U5;00->10' },
 {'name':'U5;00->10;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;00->10' , 'N':'iR:U5;00->10' , 'E':'B:U6_W_0'     },
 {'name':'U5;00->11;nw' , 'parity':1, 'S':'iL:U5;00->11' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;00->11' },
 {'name':'U5;00->11;ne' , 'parity':0, 'S':'iR:U5;00->11' , 'W':'iT:U5;00->11' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_1'     },
 {'name':'U5;00->11;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;00->11' , 'E':'iB:U5;00->11' },
 {'name':'U5;00->11;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;00->11' , 'N':'iR:U5;00->11' , 'E':'B:U6_W_1'     },
 {'name':'U5;01->00;nw' , 'parity':1, 'S':'iL:U5;01->00' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;01->00' },
 {'name':'U5;01->00;ne' , 'parity':0, 'S':'iR:U5;01->00' , 'W':'iT:U5;01->00' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_0'     },
 {'name':'U5;01->00;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;01->00' , 'E':'iB:U5;01->00' },
 {'name':'U5;01->00;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;01->00' , 'N':'iR:U5;01->00' , 'E':'B:U6_W_0'     },
 {'name':'U5;01->01;nw' , 'parity':1, 'S':'iL:U5;01->01' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;01->01' },
 {'name':'U5;01->01;ne' , 'parity':0, 'S':'iR:U5;01->01' , 'W':'iT:U5;01->01' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_1'     },
 {'name':'U5;01->01;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;01->01' , 'E':'iB:U5;01->01' },
 {'name':'U5;01->01;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;01->01' , 'N':'iR:U5;01->01' , 'E':'B:U6_W_1'     },
 {'name':'U5;01->10;nw' , 'parity':1, 'S':'iL:U5;01->10' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;01->10' },
 {'name':'U5;01->10;ne' , 'parity':0, 'S':'iR:U5;01->10' , 'W':'iT:U5;01->10' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_0'     },
 {'name':'U5;01->10;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;01->10' , 'E':'iB:U5;01->10' },
 {'name':'U5;01->10;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;01->10' , 'N':'iR:U5;01->10' , 'E':'B:U6_W_0'     },
 {'name':'U5;01->11;nw' , 'parity':1, 'S':'iL:U5;01->11' , 'W':'T:U4_E_0'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;01->11' },
 {'name':'U5;01->11;ne' , 'parity':0, 'S':'iR:U5;01->11' , 'W':'iT:U5;01->11' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_1'     },
 {'name':'U5;01->11;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_0'     , 'N':'iL:U5;01->11' , 'E':'iB:U5;01->11' },
 {'name':'U5;01->11;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;01->11' , 'N':'iR:U5;01->11' , 'E':'B:U6_W_1'     },
 {'name':'U5;10->00;nw' , 'parity':1, 'S':'iL:U5;10->00' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;10->00' },
 {'name':'U5;10->00;ne' , 'parity':0, 'S':'iR:U5;10->00' , 'W':'iT:U5;10->00' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_0'     },
 {'name':'U5;10->00;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;10->00' , 'E':'iB:U5;10->00' },
 {'name':'U5;10->00;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;10->00' , 'N':'iR:U5;10->00' , 'E':'B:U6_W_0'     },
 {'name':'U5;10->01;nw' , 'parity':1, 'S':'iL:U5;10->01' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;10->01' },
 {'name':'U5;10->01;ne' , 'parity':0, 'S':'iR:U5;10->01' , 'W':'iT:U5;10->01' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_1'     },
 {'name':'U5;10->01;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;10->01' , 'E':'iB:U5;10->01' },
 {'name':'U5;10->01;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;10->01' , 'N':'iR:U5;10->01' , 'E':'B:U6_W_1'     },
 {'name':'U5;10->10;nw' , 'parity':1, 'S':'iL:U5;10->10' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;10->10' },
 {'name':'U5;10->10;ne' , 'parity':0, 'S':'iR:U5;10->10' , 'W':'iT:U5;10->10' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_0'     },
 {'name':'U5;10->10;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;10->10' , 'E':'iB:U5;10->10' },
 {'name':'U5;10->10;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;10->10' , 'N':'iR:U5;10->10' , 'E':'B:U6_W_0'     },
 {'name':'U5;10->11;nw' , 'parity':1, 'S':'iL:U5;10->11' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_0'     , 'E':'iT:U5;10->11' },
 {'name':'U5;10->11;ne' , 'parity':0, 'S':'iR:U5;10->11' , 'W':'iT:U5;10->11' , 'N':'R:U6_S_0'     , 'E':'T:U6_W_1'     },
 {'name':'U5;10->11;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;10->11' , 'E':'iB:U5;10->11' },
 {'name':'U5;10->11;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;10->11' , 'N':'iR:U5;10->11' , 'E':'B:U6_W_1'     },
 {'name':'U5;11->00;nw' , 'parity':1, 'S':'iL:U5;11->00' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;11->00' },
 {'name':'U5;11->00;ne' , 'parity':0, 'S':'iR:U5;11->00' , 'W':'iT:U5;11->00' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_0'     },
 {'name':'U5;11->00;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;11->00' , 'E':'iB:U5;11->00' },
 {'name':'U5;11->00;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;11->00' , 'N':'iR:U5;11->00' , 'E':'B:U6_W_0'     },
 {'name':'U5;11->01;nw' , 'parity':1, 'S':'iL:U5;11->01' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;11->01' },
 {'name':'U5;11->01;ne' , 'parity':0, 'S':'iR:U5;11->01' , 'W':'iT:U5;11->01' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_1'     },
 {'name':'U5;11->01;sw' , 'parity':0, 'S':'L:U4_N_0'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;11->01' , 'E':'iB:U5;11->01' },
 {'name':'U5;11->01;se' , 'parity':1, 'S':'R:U4_N_0'     , 'W':'iB:U5;11->01' , 'N':'iR:U5;11->01' , 'E':'B:U6_W_1'     },
 {'name':'U5;11->10;nw' , 'parity':1, 'S':'iL:U5;11->10' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;11->10' },
 {'name':'U5;11->10;ne' , 'parity':0, 'S':'iR:U5;11->10' , 'W':'iT:U5;11->10' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_0'     },
 {'name':'U5;11->10;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;11->10' , 'E':'iB:U5;11->10' },
 {'name':'U5;11->10;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;11->10' , 'N':'iR:U5;11->10' , 'E':'B:U6_W_0'     },
 {'name':'U5;11->11;nw' , 'parity':1, 'S':'iL:U5;11->11' , 'W':'T:U4_E_1'     , 'N':'L:U6_S_1'     , 'E':'iT:U5;11->11' },
 {'name':'U5;11->11;ne' , 'parity':0, 'S':'iR:U5;11->11' , 'W':'iT:U5;11->11' , 'N':'R:U6_S_1'     , 'E':'T:U6_W_1'     },
 {'name':'U5;11->11;sw' , 'parity':0, 'S':'L:U4_N_1'     , 'W':'B:U4_E_1'     , 'N':'iL:U5;11->11' , 'E':'iB:U5;11->11' },
 {'name':'U5;11->11;se' , 'parity':1, 'S':'R:U4_N_1'     , 'W':'iB:U5;11->11' , 'N':'iR:U5;11->11' , 'E':'B:U6_W_1'     },
 {'name':'U3;00->00;nw' , 'parity':1, 'S':'iL:U3;00->00' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;00->00' },
 {'name':'U3;00->00;ne' , 'parity':0, 'S':'iR:U3;00->00' , 'W':'iT:U3;00->00' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_0'     },
 {'name':'U3;00->00;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;00->00' , 'E':'iB:U3;00->00' },
 {'name':'U3;00->00;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;00->00' , 'N':'iR:U3;00->00' , 'E':'B:U4_W_0'     },
 {'name':'U3;00->01;nw' , 'parity':1, 'S':'iL:U3;00->01' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;00->01' },
 {'name':'U3;00->01;ne' , 'parity':0, 'S':'iR:U3;00->01' , 'W':'iT:U3;00->01' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_1'     },
 {'name':'U3;00->01;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;00->01' , 'E':'iB:U3;00->01' },
 {'name':'U3;00->01;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;00->01' , 'N':'iR:U3;00->01' , 'E':'B:U4_W_1'     },
 {'name':'U3;00->10;nw' , 'parity':1, 'S':'iL:U3;00->10' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;00->10' },
 {'name':'U3;00->10;ne' , 'parity':0, 'S':'iR:U3;00->10' , 'W':'iT:U3;00->10' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_0'     },
 {'name':'U3;00->10;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;00->10' , 'E':'iB:U3;00->10' },
 {'name':'U3;00->10;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;00->10' , 'N':'iR:U3;00->10' , 'E':'B:U4_W_0'     },
 {'name':'U3;00->11;nw' , 'parity':1, 'S':'iL:U3;00->11' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;00->11' },
 {'name':'U3;00->11;ne' , 'parity':0, 'S':'iR:U3;00->11' , 'W':'iT:U3;00->11' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_1'     },
 {'name':'U3;00->11;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;00->11' , 'E':'iB:U3;00->11' },
 {'name':'U3;00->11;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;00->11' , 'N':'iR:U3;00->11' , 'E':'B:U4_W_1'     },
 {'name':'U3;01->00;nw' , 'parity':1, 'S':'iL:U3;01->00' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;01->00' },
 {'name':'U3;01->00;ne' , 'parity':0, 'S':'iR:U3;01->00' , 'W':'iT:U3;01->00' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_0'     },
 {'name':'U3;01->00;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;01->00' , 'E':'iB:U3;01->00' },
 {'name':'U3;01->00;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;01->00' , 'N':'iR:U3;01->00' , 'E':'B:U4_W_0'     },
 {'name':'U3;01->01;nw' , 'parity':1, 'S':'iL:U3;01->01' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;01->01' },
 {'name':'U3;01->01;ne' , 'parity':0, 'S':'iR:U3;01->01' , 'W':'iT:U3;01->01' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_1'     },
 {'name':'U3;01->01;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;01->01' , 'E':'iB:U3;01->01' },
 {'name':'U3;01->01;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;01->01' , 'N':'iR:U3;01->01' , 'E':'B:U4_W_1'     },
 {'name':'U3;01->10;nw' , 'parity':1, 'S':'iL:U3;01->10' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;01->10' },
 {'name':'U3;01->10;ne' , 'parity':0, 'S':'iR:U3;01->10' , 'W':'iT:U3;01->10' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_0'     },
 {'name':'U3;01->10;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;01->10' , 'E':'iB:U3;01->10' },
 {'name':'U3;01->10;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;01->10' , 'N':'iR:U3;01->10' , 'E':'B:U4_W_0'     },
 {'name':'U3;01->11;nw' , 'parity':1, 'S':'iL:U3;01->11' , 'W':'T:U2_E_0'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;01->11' },
 {'name':'U3;01->11;ne' , 'parity':0, 'S':'iR:U3;01->11' , 'W':'iT:U3;01->11' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_1'     },
 {'name':'U3;01->11;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_0'     , 'N':'iL:U3;01->11' , 'E':'iB:U3;01->11' },
 {'name':'U3;01->11;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;01->11' , 'N':'iR:U3;01->11' , 'E':'B:U4_W_1'     },
 {'name':'U3;10->00;nw' , 'parity':1, 'S':'iL:U3;10->00' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;10->00' },
 {'name':'U3;10->00;ne' , 'parity':0, 'S':'iR:U3;10->00' , 'W':'iT:U3;10->00' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_0'     },
 {'name':'U3;10->00;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;10->00' , 'E':'iB:U3;10->00' },
 {'name':'U3;10->00;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;10->00' , 'N':'iR:U3;10->00' , 'E':'B:U4_W_0'     },
 {'name':'U3;10->01;nw' , 'parity':1, 'S':'iL:U3;10->01' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;10->01' },
 {'name':'U3;10->01;ne' , 'parity':0, 'S':'iR:U3;10->01' , 'W':'iT:U3;10->01' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_1'     },
 {'name':'U3;10->01;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;10->01' , 'E':'iB:U3;10->01' },
 {'name':'U3;10->01;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;10->01' , 'N':'iR:U3;10->01' , 'E':'B:U4_W_1'     },
 {'name':'U3;10->10;nw' , 'parity':1, 'S':'iL:U3;10->10' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;10->10' },
 {'name':'U3;10->10;ne' , 'parity':0, 'S':'iR:U3;10->10' , 'W':'iT:U3;10->10' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_0'     },
 {'name':'U3;10->10;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;10->10' , 'E':'iB:U3;10->10' },
 {'name':'U3;10->10;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;10->10' , 'N':'iR:U3;10->10' , 'E':'B:U4_W_0'     },
 {'name':'U3;10->11;nw' , 'parity':1, 'S':'iL:U3;10->11' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_0'     , 'E':'iT:U3;10->11' },
 {'name':'U3;10->11;ne' , 'parity':0, 'S':'iR:U3;10->11' , 'W':'iT:U3;10->11' , 'N':'R:U4_S_0'     , 'E':'T:U4_W_1'     },
 {'name':'U3;10->11;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;10->11' , 'E':'iB:U3;10->11' },
 {'name':'U3;10->11;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;10->11' , 'N':'iR:U3;10->11' , 'E':'B:U4_W_1'     },
 {'name':'U3;11->00;nw' , 'parity':1, 'S':'iL:U3;11->00' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;11->00' },
 {'name':'U3;11->00;ne' , 'parity':0, 'S':'iR:U3;11->00' , 'W':'iT:U3;11->00' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_0'     },
 {'name':'U3;11->00;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;11->00' , 'E':'iB:U3;11->00' },
 {'name':'U3;11->00;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;11->00' , 'N':'iR:U3;11->00' , 'E':'B:U4_W_0'     },
 {'name':'U3;11->01;nw' , 'parity':1, 'S':'iL:U3;11->01' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;11->01' },
 {'name':'U3;11->01;ne' , 'parity':0, 'S':'iR:U3;11->01' , 'W':'iT:U3;11->01' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_1'     },
 {'name':'U3;11->01;sw' , 'parity':0, 'S':'L:U2_N_0'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;11->01' , 'E':'iB:U3;11->01' },
 {'name':'U3;11->01;se' , 'parity':1, 'S':'R:U2_N_0'     , 'W':'iB:U3;11->01' , 'N':'iR:U3;11->01' , 'E':'B:U4_W_1'     },
 {'name':'U3;11->10;nw' , 'parity':1, 'S':'iL:U3;11->10' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;11->10' },
 {'name':'U3;11->10;ne' , 'parity':0, 'S':'iR:U3;11->10' , 'W':'iT:U3;11->10' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_0'     },
 {'name':'U3;11->10;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;11->10' , 'E':'iB:U3;11->10' },
 {'name':'U3;11->10;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;11->10' , 'N':'iR:U3;11->10' , 'E':'B:U4_W_0'     },
 {'name':'U3;11->11;nw' , 'parity':1, 'S':'iL:U3;11->11' , 'W':'T:U2_E_1'     , 'N':'L:U4_S_1'     , 'E':'iT:U3;11->11' },
 {'name':'U3;11->11;ne' , 'parity':0, 'S':'iR:U3;11->11' , 'W':'iT:U3;11->11' , 'N':'R:U4_S_1'     , 'E':'T:U4_W_1'     },
 {'name':'U3;11->11;sw' , 'parity':0, 'S':'L:U2_N_1'     , 'W':'B:U2_E_1'     , 'N':'iL:U3;11->11' , 'E':'iB:U3;11->11' },
 {'name':'U3;11->11;se' , 'parity':1, 'S':'R:U2_N_1'     , 'W':'iB:U3;11->11' , 'N':'iR:U3;11->11' , 'E':'B:U4_W_1'     },
 {'name':'U8;0_->0_;nw' , 'parity':1, 'S':'iL:U8;0_->0_' , 'W':'T:U8_W_0'     , 'N':'L:U8_N'       , 'E':'iT:U8;0_->0_' },
 {'name':'U8;0_->0_;ne' , 'parity':0, 'S':'iR:U8;0_->0_' , 'W':'iT:U8;0_->0_' , 'N':'R:U8_N'       , 'E':'T:U8_E'       },
 {'name':'U8;0_->0_;sw' , 'parity':0, 'S':'L:U8_S_0'     , 'W':'B:U8_W_0'     , 'N':'iL:U8;0_->0_' , 'E':'iB:U8;0_->0_' },
 {'name':'U8;0_->0_;se' , 'parity':1, 'S':'R:U8_S_0'     , 'W':'iB:U8;0_->0_' , 'N':'iR:U8;0_->0_' , 'E':'B:U8_E'       },
 {'name':'U8;0_->1_;nw' , 'parity':1, 'S':'iL:U8;0_->1_' , 'W':'T:U8_W_0'     , 'N':'L:U8_N'       , 'E':'iT:U8;0_->1_' },
 {'name':'U8;0_->1_;ne' , 'parity':0, 'S':'iR:U8;0_->1_' , 'W':'iT:U8;0_->1_' , 'N':'R:U8_N'       , 'E':'T:U8_E'       },
 {'name':'U8;0_->1_;sw' , 'parity':0, 'S':'L:U8_S_1'     , 'W':'B:U8_W_0'     , 'N':'iL:U8;0_->1_' , 'E':'iB:U8;0_->1_' },
 {'name':'U8;0_->1_;se' , 'parity':1, 'S':'R:U8_S_1'     , 'W':'iB:U8;0_->1_' , 'N':'iR:U8;0_->1_' , 'E':'B:U8_E'       },
 {'name':'U8;1_->0_;nw' , 'parity':1, 'S':'iL:U8;1_->0_' , 'W':'T:U8_W_1'     , 'N':'L:U8_N'       , 'E':'iT:U8;1_->0_' },
 {'name':'U8;1_->0_;ne' , 'parity':0, 'S':'iR:U8;1_->0_' , 'W':'iT:U8;1_->0_' , 'N':'R:U8_N'       , 'E':'T:U8_E'       },
 {'name':'U8;1_->0_;sw' , 'parity':0, 'S':'L:U8_S_0'     , 'W':'B:U8_W_1'     , 'N':'iL:U8;1_->0_' , 'E':'iB:U8;1_->0_' },
 {'name':'U8;1_->0_;se' , 'parity':1, 'S':'R:U8_S_0'     , 'W':'iB:U8;1_->0_' , 'N':'iR:U8;1_->0_' , 'E':'B:U8_E'       },
 {'name':'U8;1_->1_;nw' , 'parity':1, 'S':'iL:U8;1_->1_' , 'W':'T:U8_W_1'     , 'N':'L:U8_N'       , 'E':'iT:U8;1_->1_' },
 {'name':'U8;1_->1_;ne' , 'parity':0, 'S':'iR:U8;1_->1_' , 'W':'iT:U8;1_->1_' , 'N':'R:U8_N'       , 'E':'T:U8_E'       },
 {'name':'U8;1_->1_;sw' , 'parity':0, 'S':'L:U8_S_1'     , 'W':'B:U8_W_1'     , 'N':'iL:U8;1_->1_' , 'E':'iB:U8;1_->1_' },
 {'name':'U8;1_->1_;se' , 'parity':1, 'S':'R:U8_S_1'     , 'W':'iB:U8;1_->1_' , 'N':'iR:U8;1_->1_' , 'E':'B:U8_E'       },
 {'name':'U2;_0->_0;nw' , 'parity':1, 'S':'iL:U2;_0->_0' , 'W':'T:U2_W'       , 'N':'L:U2_N_0'     , 'E':'iT:U2;_0->_0' },
 {'name':'U2;_0->_0;ne' , 'parity':0, 'S':'iR:U2;_0->_0' , 'W':'iT:U2;_0->_0' , 'N':'R:U2_N_0'     , 'E':'T:U2_E_0'     },
 {'name':'U2;_0->_0;sw' , 'parity':0, 'S':'L:U2_S'       , 'W':'B:U2_W'       , 'N':'iL:U2;_0->_0' , 'E':'iB:U2;_0->_0' },
 {'name':'U2;_0->_0;se' , 'parity':1, 'S':'R:U2_S'       , 'W':'iB:U2;_0->_0' , 'N':'iR:U2;_0->_0' , 'E':'B:U2_E_0'     },
 {'name':'U2;_0->_1;nw' , 'parity':1, 'S':'iL:U2;_0->_1' , 'W':'T:U2_W'       , 'N':'L:U2_N_0'     , 'E':'iT:U2;_0->_1' },
 {'name':'U2;_0->_1;ne' , 'parity':0, 'S':'iR:U2;_0->_1' , 'W':'iT:U2;_0->_1' , 'N':'R:U2_N_0'     , 'E':'T:U2_E_1'     },
 {'name':'U2;_0->_1;sw' , 'parity':0, 'S':'L:U2_S'       , 'W':'B:U2_W'       , 'N':'iL:U2;_0->_1' , 'E':'iB:U2;_0->_1' },
 {'name':'U2;_0->_1;se' , 'parity':1, 'S':'R:U2_S'       , 'W':'iB:U2;_0->_1' , 'N':'iR:U2;_0->_1' , 'E':'B:U2_E_1'     },
 {'name':'U2;_1->_0;nw' , 'parity':1, 'S':'iL:U2;_1->_0' , 'W':'T:U2_W'       , 'N':'L:U2_N_1'     , 'E':'iT:U2;_1->_0' },
 {'name':'U2;_1->_0;ne' , 'parity':0, 'S':'iR:U2;_1->_0' , 'W':'iT:U2;_1->_0' , 'N':'R:U2_N_1'     , 'E':'T:U2_E_0'     },
 {'name':'U2;_1->_0;sw' , 'parity':0, 'S':'L:U2_S'       , 'W':'B:U2_W'       , 'N':'iL:U2;_1->_0' , 'E':'iB:U2;_1->_0' },
 {'name':'U2;_1->_0;se' , 'parity':1, 'S':'R:U2_S'       , 'W':'iB:U2;_1->_0' , 'N':'iR:U2;_1->_0' , 'E':'B:U2_E_0'     },
 {'name':'U2;_1->_1;nw' , 'parity':1, 'S':'iL:U2;_1->_1' , 'W':'T:U2_W'       , 'N':'L:U2_N_1'     , 'E':'iT:U2;_1->_1' },
 {'name':'U2;_1->_1;ne' , 'parity':0, 'S':'iR:U2;_1->_1' , 'W':'iT:U2;_1->_1' , 'N':'R:U2_N_1'     , 'E':'T:U2_E_1'     },
 {'name':'U2;_1->_1;sw' , 'parity':0, 'S':'L:U2_S'       , 'W':'B:U2_W'       , 'N':'iL:U2;_1->_1' , 'E':'iB:U2;_1->_1' },
 {'name':'U2;_1->_1;se' , 'parity':1, 'S':'R:U2_S'       , 'W':'iB:U2;_1->_1' , 'N':'iR:U2;_1->_1' , 'E':'B:U2_E_1'     },
 {'name':'U1;__->__;nw' , 'parity':1, 'S':'iL:U1;__->__' , 'W':'T:U8_E'       , 'N':'L:U2_S'       , 'E':'iT:U1;__->__' },
 {'name':'U1;__->__;ne' , 'parity':0, 'S':'iR:U1;__->__' , 'W':'iT:U1;__->__' , 'N':'R:U2_S'       , 'E':'T:U2_W'       },
 {'name':'U1;__->__;sw' , 'parity':0, 'S':'L:U8_N'       , 'W':'B:U8_E'       , 'N':'iL:U1;__->__' , 'E':'iB:U1;__->__' },
 {'name':'U1;__->__;se' , 'parity':1, 'S':'R:U8_N'       , 'W':'iB:U1;__->__' , 'N':'iR:U1;__->__' , 'E':'B:U2_W'       }]



# glues that need biotins, and in which direction 
# (hence must put an internal biotinylated "T" in domain)
glue_biotins = \
[('iR:U6;00->01', 'S'),
 ('iL:U6;00->10', 'S'),
 ('iL:U6;00->11', 'S'),
 ('iR:U6;00->11', 'S'),
 ('iR:U6;01->01', 'S'),
 ('iL:U6;01->10', 'S'),
 ('iL:U6;01->11', 'S'),
 ('iR:U6;01->11', 'S'),
 ('iR:U6;10->01', 'S'),
 ('iL:U6;10->10', 'S'),
 ('iL:U6;10->11', 'S'),
 ('iR:U6;10->11', 'S'),
 ('iR:U6;11->01', 'S'),
 ('iL:U6;11->10', 'S'),
 ('iL:U6;11->11', 'S'),
 ('iR:U6;11->11', 'S'),
 ('iR:U4;00->01', 'S'),
 ('iL:U4;00->10', 'S'),
 ('iL:U4;00->11', 'S'),
 ('iR:U4;00->11', 'S'),
 ('iR:U4;01->01', 'S'),
 ('iL:U4;01->10', 'S'),
 ('iL:U4;01->11', 'S'),
 ('iR:U4;01->11', 'S'),
 ('iR:U4;10->01', 'S'),
 ('iL:U4;10->10', 'S'),
 ('iL:U4;10->11', 'S'),
 ('iR:U4;10->11', 'S'),
 ('iR:U4;11->01', 'S'),
 ('iL:U4;11->10', 'S'),
 ('iL:U4;11->11', 'S'),
 ('iR:U4;11->11', 'S'),
 ('iT:U7;00->01', 'W'),
 ('iB:U7;00->10', 'W'),
 ('iT:U7;00->11', 'W'),
 ('iB:U7;00->11', 'W'),
 ('iT:U7;01->01', 'W'),
 ('iB:U7;01->10', 'W'),
 ('iT:U7;01->11', 'W'),
 ('iB:U7;01->11', 'W'),
 ('iT:U7;10->01', 'W'),
 ('iB:U7;10->10', 'W'),
 ('iT:U7;10->11', 'W'),
 ('iB:U7;10->11', 'W'),
 ('iT:U7;11->01', 'W'),
 ('iB:U7;11->10', 'W'),
 ('iT:U7;11->11', 'W'),
 ('iB:U7;11->11', 'W'),
 ('iT:U5;00->01', 'W'),
 ('iB:U5;00->10', 'W'),
 ('iT:U5;00->11', 'W'),
 ('iB:U5;00->11', 'W'),
 ('iT:U5;01->01', 'W'),
 ('iB:U5;01->10', 'W'),
 ('iT:U5;01->11', 'W'),
 ('iB:U5;01->11', 'W'),
 ('iT:U5;10->01', 'W'),
 ('iB:U5;10->10', 'W'),
 ('iT:U5;10->11', 'W'),
 ('iB:U5;10->11', 'W'),
 ('iT:U5;11->01', 'W'),
 ('iB:U5;11->10', 'W'),
 ('iT:U5;11->11', 'W'),
 ('iB:U5;11->11', 'W'),
 ('iT:U3;00->01', 'W'),
 ('iB:U3;00->10', 'W'),
 ('iT:U3;00->11', 'W'),
 ('iB:U3;00->11', 'W'),
 ('iT:U3;01->01', 'W'),
 ('iB:U3;01->10', 'W'),
 ('iT:U3;01->11', 'W'),
 ('iB:U3;01->11', 'W'),
 ('iT:U3;10->01', 'W'),
 ('iB:U3;10->10', 'W'),
 ('iT:U3;10->11', 'W'),
 ('iB:U3;10->11', 'W'),
 ('iT:U3;11->01', 'W'),
 ('iB:U3;11->10', 'W'),
 ('iT:U3;11->11', 'W'),
 ('iB:U3;11->11', 'W'),
 ('iL:U8;0_->1_', 'S'),
 ('iL:U8;1_->1_', 'S'),
 ('iR:U2;_0->_1', 'S'),
 ('iR:U2;_1->_1', 'S')]


def direction2axis(direction):
    if direction in ['N','S']:
        return 'NS'
    if direction in ['E','W']:
        return 'EW'
    raise AssertionError('should be unreachable')

# glues that need different strengths from lowDG and highDG
glue_strength_constraints = \
            [ 
              ( 
                lowDG+biotin_strength_boost,
                highDG+biotin_strength_boost,
                [(glue_name, direction2axis(glue_direction)) for (glue_name, glue_direction) in glue_biotins]
             )
            ] 

# e.g.
                #  ("iE:U5;11->01","NS"), 
                #  ("iE:U5;11->11","NS"), 
                #  ("L:U4_N_1","NS"),


# We want to specifiy a single "zipper tile", which we do as follows (as two tiles with the same sequences):
# glues with sequences already assigned to them (note that the sequences are allowed to violate the strength constraints)
glue_sequences = [ 
#{'name':'U1_nw'   , 'parity':0, 'S':'iW:U1'      , 'W':'T:U8_E'     , 'N':'L:U2_S'     , 'E':'iN:U1'      },
                   ("iL:U1;__->__",  "S", "TGCACCTATT"),
                   ("T:U8_E", "W", "TTACTACACGA"),
                   ("L:U2_S", "N", "AAACGACAAAT"),
                   ("iT:U1;__->__",  "E", "AAGCCAATCT"),
#{'name':'U1_se'   , 'parity':0, 'S':'R:U8_N'     , 'W':'iS:U1'      , 'N':'iE:U1'      , 'E':'B:U2_W'     }]
                   ("R:U8_N", "S", "TGCACCTATT"),
                   ("iB:U1;__->__",  "W", "TTACTACACGA"),
                   ("iR:U1;__->__",  "N", "AAACGACAAAT"),
                   ("B:U2_W", "E", "AAGCCAATCT") 
                 ] 

import itertools as it



# The following variables make use of features in the sequence designer that not exploited 
# when designing DNA sequence for the 6-bit IBC tile set. In particular, assigning negative values 
# acts as a switch to cause the sequence design code to ignore the relevant parameter.

# minimum Hamming distance all domains must be from each other
hamming = -1

# threshold for allowed energy of interaction between orthogonal glues
# that will not be co-located during a strength-1 binding event
orth_any = -5.0

# orth_any is worst-case; all pairs of sequences must be below
# orth_any_ave is average orthogonality of proposed new domain with all
# currently used domains
orth_any_ave = -3.0

# this is the limit on individual secondary structure of a single domain
domain_indv_sec_struct = -0.3
domain_indv_sec_struct_ave = -0.2

# this is the limit on secondary structure of any two domains concatenated
domain_pair_sec_struct = -20.0
domain_pair_sec_struct_ave = -0.8

# If the following list is defined, it defines mutually exclusive pairs of tiles (their names, actually) that will never be together in the same pot 
# The program will not bother comparing any glues/tiles to each other if they appear as a pair here
# The code below finds pairs of tiles with the same (non-full) name (e.g., U6), 
#   and the same input to the function they compute, but different output
mutually_exclusive_tilename_pairs = []
