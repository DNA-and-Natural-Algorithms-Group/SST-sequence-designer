# parameter file to give input to sequence designer:
# uniqueNxNperiodic.py
#
# Creates an NxN array of uniquely-addressed tiles that has periodic boundariy conditions
# so that it could form an unbounded 2D array or (more likely) roll up into a tube (which would most
# likely be a 16-helix tube for N=8).
#
# The value of N is specified below; it must be even to satisfy SST tile type parity conditions.
#
# Note: The design parameters specified below are similar to those used for the complete 6-bit IBC
# tile set.  More stringent criteria could presumably be used for N=8, and presumbly they would 
# have to be relaxed for larger tile sets.

# See ibc_6bit.py for explanation of parameters defined below.


N=8
if N%2==1:
  import sys
  sys.exit("N must be even!")

lowDG=8.9
highDG=9.2
biotin_strength_boost=1.1
temperature = 53.0
tile_sec_struct = 2.0
tile_orth = 6.1
tile_orth_share = 8.1
orth_algorithmic_conflict = 1.6
orth_algorithmic_conflict_generalized = 2.6
orth_colocated = 1.4
lattice_binding_lower_threshold = 12.3
three_letter_code = True
three_letter_code_exceptions = 1
endGC = False
endAT = True
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

print("========================================================")
print("Each SST strand has domains in 5' -> 3' order:  S W N E ")
print("========================================================")
print(tiles)
print("========================================================")

glue_biotins = []

def direction2axis(direction):
    if direction in ['N','S']:
        return 'NS'
    if direction in ['E','W']:
        return 'EW'
    raise AssertionError('should be unreachable')

glue_strength_constraints = \
            [ 
              ( 
                lowDG+biotin_strength_boost,
                highDG+biotin_strength_boost,
                [(glue_name, direction2axis(glue_direction)) for (glue_name, glue_direction) in glue_biotins]
             )
            ] 
            
glue_sequences = [] 
mutually_exclusive_tilename_pairs = []

hamming = -1
orth_any = -5.0
orth_any_ave = -3.0
domain_indv_sec_struct = -0.3
domain_indv_sec_struct_ave = -0.2
domain_pair_sec_struct = -20.0
domain_pair_sec_struct_ave = -0.8
domain_pair_sec_struct_ave = -0.8
