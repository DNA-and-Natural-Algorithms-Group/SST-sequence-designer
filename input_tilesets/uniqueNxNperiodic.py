# parameter file to give input to sequence designer:
# uniqueNxNperiodic.py
#
# Creates an NxN array of uniquely-addressed tiles that has periodic boundariy conditions
# so that it could form an unbounded 2D array or (more likely) roll up into a tube.
#
# The value of N is specified below; it must be even to satisfy SST tile type parity conditions.

N=6
if N%2==1:
  import sys
  sys.exit("N must be even!")


# Note that all energy parameter values are sign-flipped 
# (i.e. more positive is more favourable)


# "global" energy interval for all glues, except those glues with user-specified sequences
# evaluated using a simple nearest-neighbour model
lowDG=8.9
highDG=9.2
biotin_strength_boost=1.1

# temperature in degrees C at which to evaluate energies
temperature = 53.0

# allowed energy of secondary structure of each individual tile (via nupack)
tile_sec_struct = 2.0

# allowed energy of secondary structure of each pair of tiles that don't share any glues (via RNAduplex)
tile_orth = 6.1

# allowed energy of secondary structure of pair of tiles that share glues (via RNAduplex)
tile_orth_share = 8.1

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






# algorithmically build a NxN set of tiles, each tile is a python dictionary
tiles=list()
for i in range(N):
  for j in range(N):
    tiles.append({
            'name': 'T'+str(i)+'x'+str(j),
            'parity': (i+j)%2,
            'N': 'glue_NS_'+str((i+1)%N)+'x'+str(j),
            'E': 'glue_EW_'+str(i)+'x'+str((j+1)%N),
            'S': 'glue_NS_'+str(i)+'x'+str(j),
            'W': 'glue_EW_'+str(i)+'x'+str(j)
    })

print "========================================================"
print "Each SST strand has domains in 5' -> 3' order:  S W N E "
print "========================================================"
print tiles
print "========================================================"


# glues that need biotins, and in which direction 
# (hence must put an internal biotinylated "T" in domain)
glue_biotins = []

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
            

# Glues with sequences already assigned to them (note that the sequences are allowed to violate the strength constraints).
# In particular this is useful when specifiying a single "seam tile" for two posiitons in a proof-reading block, 
# which we do here as follows (as two tiles with the same sequences):
glue_sequences = [] 

# If this is defined, it defines mutually exclusive pairs of tiles (their names, actually) 
# that will never be together in the same pot 
# The program will not bother comparing any glues/tiles to each other if they appear as a pair here
# The code below finds pairs of tiles with the same (non-full) name (e.g., U6), 
#   and the same input to the function they compute, but different output
mutually_exclusive_tilename_pairs = []







# Below, negative values cause code to ignore the relevant parameter
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
domain_pair_sec_struct_ave = -0.8
