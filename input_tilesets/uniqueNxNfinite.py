# parameter file to give input to sequence designer:
# uniqueNxNfinite.py
#
# Creates an NxN array of uniquely-addressed tiles that has poly-T domains on the boundary.
# 
# Note 1: While the sequence designer will produce tile sequences for this tile set, please 
# be aware that atam2ssts was formulated with DNA nanotubes in mind, and with a specific
# growth direction, and thus we have not evaluated the suitability of these sequences for
# other growth patterns.  In fact, we recommend that they be taken with extreme caution.
#
# Note 2: Because the poly-T sequences are assigned to several glues with different names,
# the "lattice-binding" score will never reach zero and thus the designer will never halt.
# You will have to use control-C to stop the process, and then manually kill any residual
# pfunc processes that may be left running.  But the sequences left in the output file,
# which is written whenever a new "best sequence" is found, should be usable.  
#
# Note 3: You can avoid using the poly-T sequences by setting "glue_sequences = []" and
# in that case the design process should come to a natural halt.
#
# Note 4: This NxN finite array has a similar layout to to the 310-pixel rectangular canvas of 
# "Complex shapes self-assembled from single-stranded DNA tiles" (Wei, Dai, Yin, Nature, 2012)
# but has different boundary conditions and orientation (diamond vs square).

# See ibc_6bit.py for explanation of parameters defined below.

N=9

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
            'N': 'glue_NS_'+str(i+1)+'x'+str(j),
            'E': 'glue_EW_'+str(i)+'x'+str(j+1),
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
            
glue_sequences = [ ("glue_NS_%dx%d"%(0,i),"S","TTTTTTTTTTT" if i%2==0 else "TTTTTTTTTT") for i in range(N) ] + \
                 [ ("glue_NS_%dx%d"%(N,i),"N","TTTTTTTTTTT" if (N+i)%2==0 else "TTTTTTTTTT") for i in range(N) ] + \
                 [ ("glue_EW_%dx%d"%(i,0),"W","TTTTTTTTTTT" if i%2==1 else "TTTTTTTTTT") for i in range(N) ] + \
                 [ ("glue_EW_%dx%d"%(i,N),"E","TTTTTTTTTTT" if (N+i)%2==1 else "TTTTTTTTTT") for i in range(N) ]

print glue_sequences

mutually_exclusive_tilename_pairs = []



hamming = -1
orth_any = -5.0
orth_any_ave = -3.0
domain_indv_sec_struct = -0.3
domain_indv_sec_struct_ave = -0.2
domain_pair_sec_struct = -20.0
domain_pair_sec_struct_ave = -0.8
domain_pair_sec_struct_ave = -0.8
