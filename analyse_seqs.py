'''Python 3 script to analyse DNA single-stranded tile sequences.'''


#import matplotlib.patches as mpatches
#import multiprocessing
# string,random,


import math,sys,os,time,itertools,pickle,pprint  
import subprocess as sub
import numpy as np
import argparse as argparse

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing.pool import ThreadPool
from functools import lru_cache 
from datetime import datetime
global_thread_pool = ThreadPool()


MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE = True 


viennaRNA_PARAMETER_SET_DIRECTORY = 'nupack_viennaRNA/'
DEFAULT_viennaRNA_PARAMETER_SET = viennaRNA_PARAMETER_SET_DIRECTORY+'dna_mathews1999.par'


#############################################################################
# Define constants and inputs  
#############################################################################
_default_temp_in_C = 53.0 # if the user does not specify a temperature, we use this value
# _idt_file = True
_threaded = True 
_bad_first_algo_conflict_threshold = -8.0
_lower_thres_correct_binding = -10.8

_check_double_mismatches = False # Setting to False speeds up lattice binding analysis
_num_helices=16    
_bad_lattice_binding_to_arbitrary_lattice_threshold = -10.0

# Some output filenames
filelist = ['nupack domain pair binding V RNAduplex domain pair (1999) - non-wc.pdf',
            'nupack domain pair binding V RNAduplex domain pair (1999) - wc.pdf',
            'dG of lattice+tile-lattice-tile V binding(tile,lattice).pdf',
            'nupack pair binding V RNAduplex pair (1999).pdf',
            'nupack pair binding V RNAduplex pair (2004).pdf',
            'RNAduplex pair (1999) V RNAduplex pair (2004).pdf',
            'nupack pair mfe V RNAduplex pair.pdf',
            'nupack pair binding V nupack pair mfe.pdf']


#############################################################################
# args e.g. --help  
#############################################################################

def parse_args():
  parser = argparse.ArgumentParser(description='''Python3 script that runs a barrage of energetics 
  tests on DNA-tile sequences. If you are running a non-algorithmic tile set please use option -nl 
  (to turn off some tests and prevent errors). Input is given as a text file in IDT format, for example:
  
  # Example 1 
  # Algorithmic SST example. Note that sequences may contain T bases with an ("internal") biotin modification, denoted "/iBiodT/". 
  # Strand domains (that encode aTAM tile-glues) are separated by a single space. Comments are ignored by the analysis code.
  # 
  U3;10->10;sw,AGTGTGTTTTT AGTTCGATGT AGAGGCTTTT ATCAGAGGGAA,25nm,STD
  U3;10->10;se,ACTCGTTCTT TTCCC/iBiodT/CTGAT TCCTCCAATTA AAAAACGCAA,25nm,STD
  U3;10->11;nw,AAGAACCACT TTCGTCAAATT ATACATCACCT AACCCACCAA,25nm,STD

  # Example 2
  # Non-algorithmic example, tile names do not have any ";" symbols. 
  # This example shows that tiles need not have four domains/glues. 
  # 
  tile_name1,AATCCTAGAA ATTGTTATTTC ATGTATACAAA GATAGATCAG,25nm,STD
  tile_name2,AACTAGAAAC CTTAGGAATT,25nm,STD
  tile_name3,TAATACTTTCA TTTATCATCG CTACATTCTT ATTTGTTTATC,25nm,STD

  The script makes use of pfunc and mfe from NUPACK, and RNAduplex from ViennaRNA, and assumes these
  executables are available. 
''',
formatter_class=argparse.RawDescriptionHelpFormatter,
epilog= '''Example usage:
python analyse_seqs.py input_file.idt
''')
  #parser.add_argument('-s', '--seqs', required=True, # action="store_true",
  #        help="Specify a batch file that contains one or more image filenames to process.", default="")
  # parser.add_argument('-o', '--check', required=False, action="store_true",  # action="store", dest="query",
  #         help="Use with --batch. Check orientation of ribbons. Produces ribbon image files, some of which may be incorrectly orientated.", default="")

  parser.add_argument('files', nargs='*') # default=os.getcwd()
  #_StoreAction(option_strings=[], dest='files', nargs='*', const=None, default='', type=None, choices=None, help=None, metavar=None)

  #parser.add_argument('-s', '--strands', required=False, action="store_true",  # action="store", dest="query",
  #        help="Analyse secondary structure energy of individual strands.", default="")

  parser.add_argument('-T', '--temperature', required=False, type=float,  # action="store", dest="query",
          help="Temperature in C to analyse sequences. Default temperature is 53.0 C")
  parser.add_argument('-nl', '--no-lattice-binding', required=False, action="store_true",  # action="store", dest="query",
          help="Do not run analysis of lattice-binding properties (includes analysis of algorithmic errors -- probably only works for algorithmic tile sets.", default="")
  parser.add_argument('-nd', '--no-domains', required=False, action="store_true",  # action="store", dest="query",
          help="Do not run domain-level secondary structure analysis.", default="")
  parser.add_argument('-ns', '--no-strand-analysis', required=False, action="store_true",  # action="store", dest="query",
          help="Do not run strand-level/tile-level secondary structure analysis.", default="")
  parser.add_argument('-np', '--no-tile-pairs', required=False, action="store_true",  # action="store", dest="query",
          help="Do not run the (expensive) tile pair test, nor the domain pair test.", default="")
  parser.add_argument('-nm', '--no-mfe-analysis', required=False, action="store_true",  # action="store", dest="query",
          help="Do not run (expensive) mfe analysis.", default="")
  parser.add_argument('-p', '--pickle', required=False, action="store_true",  # action="store", dest="query",
          help="Use python pickle lbrary to save processed data. Useful when doing time-consuming calculations and want the data for further processing.", default="")
  parser.add_argument('-b', '--bins', required=False, type=int,  # action="store", dest="query",
          help="Number of bins in histogram plots.")

  parser.add_argument('-nst', '--non-standard-tiles', required=False, action="store_true",  # action="store", dest="query",
          help="Tiles that do not have exactly 4 domains.", default="")

  args = parser.parse_args()
  if not (args.files):
    parser.error('Error: no input filename specified. Please specify a .idt file with DNA sequences for processing.')
  return args



#############################################################################
# Run analysis  
#############################################################################

#_wctable = string.maketrans('ACGTacgt_-','TGCAtgca-_')
#def wc(seq):
#    '''Return reverse Watson-Crick complement of seq'''
#    return seq.translate(_wctable)[::-1]

_wctable = bytes.maketrans(b'ACGTacgt_-',b'TGCAtgca-_')
def wc(seq):
    '''Return reverse Watson-Crick complement of seq'''
    return seq.translate(_wctable)[::-1]


def run_analysis(f,directory,temp_in_C,domain_analysis=1,strand_analysis=1,tile_pair_check=1,
                lattice_binding_analysis=1,run_mfe=1,non_standard_tiles=0,pickle_data=0,num_bins=-1):
  
  domains, domains_with_biotin = read_domains_from_idt_order(f,non_standard_tiles)
  names, seqs, seqs_bt, seq_ws, seq_ws_bt = read_strand_data_from_idt_order(f)

  if not tile_pair_check: num_tile_and_domain_pairs = -1  # -1 means do not do the tile pair, nor domain pair, check
  else: num_tile_and_domain_pairs = 0 # 0 means do the tile pair, and domain pair, check (for fast debugging/testing: num_tile_and_domain_pairs>1 means only check the first num_tile_and_domain_pairs pairs)
  
  # strand level analysis
  if strand_analysis:
    run_strand_analysis(seqs,seqs_bt,domains,domains_with_biotin,temp_in_C,
              num_tile_pairs=num_tile_and_domain_pairs, run_mfe=mfe_analysis,
              directory=directory,pickle_data=pickle_data,num_bins=num_bins)

  if domain_analysis:
    # domain level analysis
    analyse_domains(names,seqs,seqs_bt,domains,domains_with_biotin,
              seq_ws,seq_ws_bt,temp_in_C,num_domain_pairs=num_tile_and_domain_pairs,
              directory=directory,non_standard_tiles=non_standard_tiles,pickle_data=pickle_data,num_bins=num_bins) # num_tile_pairs  ,num_domain_pairs=num_tile_pairs


  if lattice_binding_analysis and not non_standard_tiles:
    print('\nStarting lattice binding analysis')
    list_of_pairs_of_lists,descriptors = run_lattice_binding_analysis(f,directory,temp_in_C)  
    plot_energies2(list_of_pairs_of_lists, labels=descriptors,
                   xaxis='dG of lattice+tile-lattice-tile',yaxis='binding(tile,lattice)', 
                   title_info='',plot_dir=directory)


  #convert_scatter_plots_to_jpg(f, directory) 
  #print('Done with: run_analysis()')


# analyze sequences at the 'strand' level of abstraction
def run_strand_analysis(seqs,seqs_bt,domains,domains_with_biotin,temp_in_C,
  num_tile_pairs=False,run_mfe=False,directory='plots/',pickle_data=0,num_bins=-1):  
    '''Analyse `strand' energetics: secondary structure, and strand pair interactions'''
    print('running analysis on strands')
    description = '' # f[1]
    
    if not os.path.exists(directory):
        print(directory)
        os.makedirs(directory)

    pickle_directory = directory+'pickle/'
    if pickle_data:
      if not os.path.exists(pickle_directory):
        os.makedirs(pickle_directory)
    
    print('  run_strand_analysis(): plots will be placed in: '+directory+  pickle_data*(', pickled datastructures will be placed in:  ' + pickle_directory))  

    print(str(len(seqs_bt))+' sequences (where we leave in biotin markers), ' + str(len(seqs))+' sequences (where we do not leave in biotin markers)')
    print(str(len([s for s in seqs_bt if '/iBiodT/' in  s])) +' sequences with biotin markers')
    if all([len(s)==42 for s in seqs] ): print('All strands are of length 42 bases')
    else: 'Some strands are not of length 42 bases'
    print(str(len(domains))+' unique domains (glue sequences), ' +str(len([d for d in domains if len(d)==10])) + ' of length 10, '+str(len([d for d in domains if len(d)==11])) + ' of length 11')
    sys.stdout.flush()

    # domain secondary structure
    #print 'analysing domain secondary structure';sys.stdout.flush()
    #print temp_in_C
    #domain_energies = [pfunc(d, temp_in_C) for d in domains] 
    #histogram([domain_energies],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='domain secondary structure',plot_dir=directory)
    #pickle.dump(domain_energies, open( pickle_directory+"domain_secondary_structure.p", "wb" ) )
    
    # strand secondary structure
    sys.stdout.write('analysing tile (strand) secondary structure: '); sys.stdout.flush()
    sec_struct_energies = [pfunc(s, temp_in_C) for s in seqs] 
    sys.stdout.write(str(len(sec_struct_energies))+' tiles\n'); sys.stdout.flush()
    histogram([sec_struct_energies],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='internal tile secondary structure (via Nupack pfunc)',plot_dir=directory,num_bins=num_bins)
    
    if pickle_data:
      pickle.dump(sec_struct_energies, open(pickle_directory+"sec_struct_energies.p", "wb" ) )
    
    # tile pair analysis
    if num_tile_pairs != -1:
      sys.stdout.write('analysing tile pairs: ');sys.stdout.flush()
      seq_pairs = get_seq_pairs(seqs)
      if num_tile_pairs>0 : seq_pairs = seq_pairs[:num_tile_pairs]
      print(str(len(seq_pairs)) +  ' tile pairs to process');sys.stdout.flush()
    
      RNAduplex_pair_energies = tile_pair_RNAduplex_energies(seq_pairs,temp_in_C)
      RNAduplex_pair_energies_2004 = tile_pair_RNAduplex_energies(seq_pairs,temp_in_C,NA_parameter_set='dna_mathews2004.par')
      histogram([RNAduplex_pair_energies],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='tile pairs - RNAduplex (1999)',plot_dir=directory,num_bins=num_bins)
      histogram([RNAduplex_pair_energies_2004],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='tile pairs - RNAduplex (2004)',plot_dir=directory,num_bins=num_bins)
    
      binding_pair_energies = tile_pair_NUPACK_binding_energies(seq_pairs,temp_in_C)
      histogram([binding_pair_energies],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='tile pairs - nupack binding',plot_dir=directory,num_bins=num_bins)
      scatter_plot_energies(binding_pair_energies, RNAduplex_pair_energies, xaxis='nupack pair binding', yaxis='RNAduplex pair (1999)', title_info='',plot_dir=directory) #,xthreshold=9.3,ythreshold=9.3
      scatter_plot_energies(binding_pair_energies, RNAduplex_pair_energies_2004, xaxis='nupack pair binding', yaxis='RNAduplex pair (2004)', title_info='',plot_dir=directory) #,xthreshold=9.3,ythreshold=9.3
      scatter_plot_energies(RNAduplex_pair_energies, RNAduplex_pair_energies_2004, xaxis='RNAduplex pair (1999)', yaxis='RNAduplex pair (2004)', title_info='',plot_dir=directory)
        
    if run_mfe and num_tile_pairs != -1:
      mfe_pair_energies = tile_pair_NUPACK_mfe_energies(seq_pairs,temp_in_C)
      histogram([mfe_pair_energies],labels=[description],xaxis='energy (kcal/mol)',yaxis='',title='tile pairs - nupack mfe',plot_dir=directory,num_bins=num_bins)
      scatter_plot_energies(mfe_pair_energies, RNAduplex_pair_energies, xaxis='nupack pair mfe', yaxis='RNAduplex pair', title_info='',plot_dir=directory)
      scatter_plot_energies(binding_pair_energies, mfe_pair_energies, xaxis='nupack pair binding', yaxis='nupack pair mfe', title_info='',plot_dir=directory)
      descriptors = ['RNAduplex (1999)','RNAduplex (2004)','NUPACK binding','NUPACK mfe']
      histogram([RNAduplex_pair_energies,RNAduplex_pair_energies_2004,binding_pair_energies,mfe_pair_energies],labels=descriptors,xaxis='energy (kcal/mol)',yaxis='',title='tile pair energies',plot_dir=directory,num_bins=num_bins)
      if pickle_data:
        pickle.dump(mfe_pair_energies, open( pickle_directory+"mfe_pair_energies.p", "wb" ) )
    elif not run_mfe and num_tile_pairs != -1:
      descriptors = ['RNAduplex (1999)','RNAduplex (2004)','NUPACK binding']
      histogram([RNAduplex_pair_energies,RNAduplex_pair_energies_2004,binding_pair_energies],labels=descriptors,xaxis='energy (kcal/mol)',yaxis='',title='tile pair energies',plot_dir=directory,num_bins=num_bins)
        
    #descriptors = ['binding v RNAduplex (1999)','binding v RNAduplex (2004)','mfe v RNAduplex (1999)','mfe v RNAduplex (2004)']
    
    #pickle.dump(RNAduplex_pair_energies, open( pickle_directory+"RNAduplex_pair_energies.p", "wb" ) )
    #pickle.dump(binding_pair_energies, open( pickle_directory+"binding_pair_energies.p", "wb" ) )
    #pickle.dump(seqs, open( pickle_directory+"seqs.p", "wb" ) )
    #pickle.dump(seq_pairs, open( pickle_directory+"seq_pairs.p", "wb" ) )
    
    #if run_mfe: return seqs, seq_pairs, sec_struct_energies, binding_pair_energies, RNAduplex_pair_energies, binding_pair_energies, mfe_pair_energies
    #else:       return seqs, seq_pairs, sec_struct_energies, binding_pair_energies, RNAduplex_pair_energies, binding_pair_energies



def lattice_binding_spacer(strand,lattice_end_top,lattice_spacer,lattice_end_bottom,temp_in_C):
    '''Strand is string with 4 single-whitespace delimited domains. lattice is of the form domain+lattice_spacer+domain,
    where lattice_spacer is typically 'TTTTT'.'''
    tile1_domains = strand.split(' ')
    dG_sticky_ends_bound = max(binding(tile1_domains[1]+tile1_domains[2],lattice_end_top+   lattice_spacer +lattice_end_bottom,temp_in_C),
                               binding(tile1_domains[1]+tile1_domains[2],lattice_end_top+wc(lattice_spacer)+lattice_end_bottom,temp_in_C))   
    dG_tube_ends_after   = max(pfunc(tile1_domains[0]+   lattice_spacer +tile1_domains[3],temp_in_C),
                               pfunc(tile1_domains[0]+wc(lattice_spacer)+tile1_domains[3],temp_in_C))
    dG_strand_before     = pfunc(''.join(tile1_domains),temp_in_C)
    dG_tube_ends_before  = max(pfunc(lattice_end_top+   lattice_spacer +lattice_end_bottom,temp_in_C),
                               pfunc(lattice_end_top+wc(lattice_spacer)+lattice_end_bottom,temp_in_C))
    dG = dG_sticky_ends_bound + dG_tube_ends_after - dG_strand_before - dG_tube_ends_before # + dG_adjustment
    return dG

_bad_single_mismatch_threshold = 10.0


def run_lattice_binding_analysis(f,directory,temp_in_C,lattice_spacer='TTTTT',num_bins=-1):    
  # check for "first" algo conflicts (single mismatches that are the first error in a proofreading block)
  # first with empty lattice_spacer
  _,_,_ = analyse_algo_conflicts(f,directory,temp_in_C=temp_in_C,lattice_spacer='',threaded=True,num_bins=num_bins) 

  # then with lattice_spacer='TTTTT'
  _,_,_ = analyse_algo_conflicts(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=True,num_bins=num_bins) 
  _,_,_ = analyse_algo_conflicts(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=False,num_bins=num_bins) 
  correct_dG_algo,mismatch_below_algo,mismatch_above_algo = analyse_algo_conflicts(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=True) 
  
  # check for generalised algo conflicts (single mismatches)
  analyse_row_conflicts(f,directory, temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,normed=True,check_double_mismatches=_check_double_mismatches,num_bins=num_bins)
  analyse_row_conflicts(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,normed=False,check_double_mismatches=_check_double_mismatches,num_bins=num_bins)
  
  _,_,_ = lattice_binding_energies_detailed(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=True,num_bins=num_bins)
  _,_,_ = lattice_binding_energies_detailed(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=False,num_bins=num_bins)
  correct_dG,mismatch_below_dG,mismatch_above_dG = lattice_binding_energies_detailed(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,threaded=True,num_bins=num_bins)
  
  correct_binding,mismatch_below_binding,mismatch_above_binding = lattice_binding_energies_simplified(f,directory,temp_in_C=temp_in_C,lattice_spacer=lattice_spacer,num_bins=num_bins)

  list_of_pairs_of_lists = [[correct_dG,correct_binding],
                            [mismatch_below_dG,mismatch_below_binding],
                            [mismatch_above_dG,mismatch_above_binding]]
  descriptors = ['correct binding',
                  '1 mismatch (below)',
                  '1 mismatch (above)']
  return list_of_pairs_of_lists,descriptors


def algo_format_tile_name(n):
  return ( len(n.split(";"))==3  and  n.split(";")[2] in ['ne','nw','se','sw'])


def analyse_row_conflicts(f,directory,temp_in_C,lattice_spacer='TTTTT',check_double_mismatches=False,normed=True,threaded=True,num_bins=-1):
  '''Nupack binding() energy of row-conflicting tile (single-mismatch that can arise 
  during correct and incorrect growth) ... TODO: explain row-conflict'''
  
  names, _, _, seqs, seqs_bt  = read_strand_data_from_idt_order(f)
  domains, domains_incl_biotin_labels = read_domains_from_idt_order(f)


  is_algo_tile_set = True
  for n in names:
    if not algo_format_tile_name(n): 
      is_algo_tile_set = False
  if not(is_algo_tile_set): 
    print("The input tile set seems to not be an algorithmic tile set (based on it's tile names), \
hence I'm skipping execution of analyse_row_conflicts().")
    return 0

  if not os.path.exists(directory): os.makedirs(directory)
  print(' analyse_row_conflicts(), lattice binding analysis from file:  '+f+'\n plots will be placed in: ' + directory)

  names_seqs = list(zip(names,seqs))

  #threaded = False # _threaded # threaded code is not yet implemented
  num_tiles_processed = 0
  correct_growth_energies = []
  single_mismatch_below_energies = []
  single_mismatch_above_energies = []
  double_mismatch_energies = []
  double_mismatch_crazy_energies = []
  thres = _bad_single_mismatch_threshold

  helix_tiles = [[] for i in range(_num_helices+1)]

  for helix in range(_num_helices):
    for tile_name, tile_seq in names_seqs:
      if tile_helix(tile_name) == helix:
        helix_tiles[helix].append( (tile_name, tile_seq) )

  for tile1_name, tile1_seq in names_seqs: 
    if num_tiles_processed%10==0 and num_tiles_processed>0: sys.stdout.write(str(num_tiles_processed))
    num_tiles_processed += 1
    sys.stdout.write('.'); sys.stdout.flush()
  
    correct_growth_energies.append(lattice_binding_spacer(tile1_seq,
                                                       wc(tile1_seq.split(' ')[2]),
                                                       lattice_spacer,
                                                       wc(tile1_seq.split(' ')[1]),
                                                       temp_in_C))
      
    # mismatch above:
    for tile2_name,tile2_seq in helix_tiles[ (tile_helix(tile1_name)+1) % _num_helices ]:
      # i.e. both tile1 and tile2 are "e"ast or both are "w"est within their pf blocks
      if tile_dir(tile1_name)[1] == tile_dir(tile2_name)[1] and  tile1_seq.split(' ')[2] != wc(tile2_seq.split(' ')[0]):  # really I should check tile names intead of seqs here                                       
        e = lattice_binding_spacer(tile1_seq,
                                  tile2_seq.split(' ')[0],
                                  lattice_spacer,
                                  wc(tile1_seq.split(' ')[1]),
                                  temp_in_C)
        single_mismatch_above_energies.append(e)
        if e > thres: 
          print_lattice_mismatch_above(e, tile2_seq.split(' ')[0], tile2_name, tile1_seq, tile1_name, wc(tile1_seq.split(' ')[1]) )
      
    # mismatch below:
    h = (tile_helix(tile1_name)-1) % _num_helices 
    if h == 0: h = 16
      
    for tile2_name,tile2_seq in helix_tiles[ h ]:
      # i.e. one of tile1, tile2 is "e"ast and the other is "w"est within their pf blocks
      if tile_dir(tile1_name)[1] != tile_dir(tile2_name)[1] and tile1_seq.split(' ')[1] != wc(tile2_seq.split(' ')[3]): # really I should check tile names intead of seqs here         
        e = lattice_binding_spacer(tile1_seq,
                                  wc(tile1_seq.split(' ')[2]),
                                  lattice_spacer,
                                  tile2_seq.split(' ')[3],
                                  temp_in_C)
        single_mismatch_below_energies.append(e)
        if e > thres: 
          print_lattice_mismatch_below(e, tile1_seq.split(' ')[2], tile1_seq, tile1_name, wc(tile2_seq.split(' ')[3]), tile2_name)
    list_of_lists_of_energies= [correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies]
    labels=['correct growth ({} energies)'.format(len(correct_growth_energies)),
            'one mismatch below ({} energies)'.format(len(single_mismatch_below_energies)),
            'one mismatch above ({} energies)'.format(len(single_mismatch_above_energies))]
      
    # Mismatch below and above (double mismatch)
    if check_double_mismatches:
      # double mismatches that do preserve row/direction
      for tile_below_name,tile_below_seq in helix_tiles[ h ]: # h was defined above             
        double_mismatch_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
          tile_above_seq.split(' ')[0],
          lattice_spacer,
          tile_below_seq.split(' ')[3],
          temp_in_C)) for tile_above_name,tile_above_seq in helix_tiles[ (tile_helix(tile1_name)+1) % _num_helices ] \
          if tile1_seq.split(' ')[1] != wc(tile_below_seq.split(' ')[3]) and \
            tile1_seq.split(' ')[2] != wc(tile_above_seq.split(' ')[0]) and \
            tile_dir(tile1_name)[1] != tile_dir(tile_below_name)[1] and \
            tile_dir(tile1_name)[1] == tile_dir(tile_above_name)[1]]        
        double_mismatch_energies.extend([result.get() for result in double_mismatch_res]) 

      # double mismatches that do not preserve row/direction
      for tile_below_name,tile_below_seq in helix_tiles[ h ]: # h was defined above             
        double_mismatch_crazy_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
            tile_above_seq.split(' ')[0],
            lattice_spacer,
            tile_below_seq.split(' ')[3],
            temp_in_C)) for tile_above_name,tile_above_seq in helix_tiles[ (tile_helix(tile1_name)+1) % _num_helices ] \
            if tile1_seq.split(' ')[1] != wc(tile_below_seq.split(' ')[3]) and \
                tile1_seq.split(' ')[2] != wc(tile_above_seq.split(' ')[0])]        
        double_mismatch_crazy_energies.extend([result.get() for result in double_mismatch_crazy_res]) 

  if check_double_mismatches:
    list_of_lists_of_energies.extend( [double_mismatch_crazy_energies, double_mismatch_energies] )
    labels.extend(['double mismatch, arbitrary ({} energies)'.format(len(double_mismatch_crazy_energies)),
                  'double mismatch, row/dir preserving ({} energies)'.format(len(double_mismatch_energies))])
                 
  histogram(list_of_lists_of_energies,labels=labels,
            xaxis='energy (kcal/mol)',yaxis='',
            title='growth on right (dG before-after binding, 1'+ ' and 2'*check_double_mismatches +' mismatch'+'es'*check_double_mismatches+', lattice spacer is '+lattice_spacer*(bool(lattice_spacer))+'empty string'*(not(bool(lattice_spacer)))+', not normalised'*(not normed)+')', #; for {} tiles)'.format(len(names_seqs)),
            label_size=12,legend_font_size=9,
            plot_dir=directory,normed=normed,xmin=-18,num_bins=num_bins) # ,xmin=4,xmax=18
  return correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies,double_mismatch_energies


def lattice_binding_energies_detailed(f,directory,temp_in_C=53,lattice_spacer='TTTTT',threaded=True,num_bins=-1):
  '''Nupack binding() energy of tile binding of right (non-seeded end) of 
  nanotube with zero or one mistmatches'''

  names, _, _, seqs, seqs_bt  = read_strand_data_from_idt_order(f)
  domains, domains_incl_biotin_labels = read_domains_from_idt_order(f)
  if not os.path.exists(directory): os.makedirs(directory)
  print('lattice_binding_energies_detailed(), lattice binding analysis from file:\n'+f+' plots will be placed in: ' + directory)
  names_seqs = list(zip(names,seqs))

  num_tiles_processed = 0
  correct_growth_energies = []
  single_mismatch_below_energies = []
  single_mismatch_above_energies = []
  bad_conflicts = ''

  for tile1_name, tile1_seq in names_seqs:   
    tile1_domains = tile1_seq.split(' ')
    if threaded:
      if num_tiles_processed%10==0 and num_tiles_processed>0: sys.stdout.write(str(num_tiles_processed))
      sys.stdout.write('.'); sys.stdout.flush()
      num_tiles_processed += 1

      correct_growth_energies.append(lattice_binding_spacer(tile1_seq,
                                                     wc(tile1_domains[2]),
                                                     lattice_spacer,
                                                     wc(tile1_domains[1]),
                                                     temp_in_C))
  
      mismatch_below_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
                                                           #wc(tile1_domains[2])+'TTTTT'+tile2_seq.split(' ')[3],
                                                           wc(tile1_domains[2]),
                                                           lattice_spacer,
                                                           tile2_seq.split(' ')[3],
                                                           temp_in_C))  for _,tile2_seq in names_seqs if wc(tile1_seq.split(' ')[1])!=tile2_seq.split(' ')[3]]
  
      single_mismatch_below_energies = single_mismatch_below_energies + [result.get() for result in mismatch_below_res] 

      mismatch_above_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
                                                           #tile2_seq.split(' ')[0]+'TTTTT'+wc(tile1_domains[1]),
                                                           tile2_seq.split(' ')[0],
                                                           lattice_spacer,
                                                           wc(tile1_domains[1]),
                                                           temp_in_C)) for _,tile2_seq in names_seqs if wc(tile1_seq.split(' ')[2])!=tile2_seq.split(' ')[0] ]
      single_mismatch_above_energies = single_mismatch_above_energies + [result.get() for result in mismatch_above_res] 

    else: # not threaded
      thres = _bad_lattice_binding_to_arbitrary_lattice_threshold 
      correct_growth_energies.append(lattice_binding_spacer(tile1_seq,
                                                     wc(tile1_domains[2]),
                                                     lattice_spacer,
                                                     wc(tile1_domains[1]),
                                                     temp_in_C))
      
      for tile2_name,tile2_seq in names_seqs:
        # mismatch below
        if wc(tile1_domains[1])!=tile2_seq.split(' ')[3]:
          e = lattice_binding_spacer(tile1_seq,
                              wc(tile1_domains[2]),
                              lattice_spacer,
                              tile2_seq.split(' ')[3],
                              temp_in_C)
          single_mismatch_below_energies.append(e)
          if e < thres: 
            bad_conflicts += print_lattice_mismatch_below(e, wc(tile1_seq.split(' ')[2]), 
                                                          tile1_seq, tile1_name, 
                                                          tile2_seq.split(' ')[3], 
                                                          'mismatch: dom 3 of '+tile2_name)
            bad_conflicts += str('{0:.2f}'.format(eval_colocated_end_pair(tile1_seq.split(' ')[1],tile2_seq.split(' ')[3],temp_in_C))) +' = eval_colocated_end_pair('+str( tile1_seq.split(' ')[1])+','+str(tile2_seq.split(' ')[3])+','+str(temp_in_C)+')\n'
                
          # mismatch above
          if wc(tile1_domains[2])!=tile2_seq.split(' ')[0]:
            e = lattice_binding_spacer(tile1_seq,
                                tile2_seq.split(' ')[0],
                                lattice_spacer,
                                wc(tile1_domains[1]),
                                temp_in_C)
            single_mismatch_above_energies.append(e)
            if e < thres:
              bad_conflicts += print_lattice_mismatch_above(e, tile2_seq.split(' ')[0], 
                                                            'mismatch: dom 0 of tile '+tile2_name, 
                                                            tile1_seq, tile1_name, 
                                                            wc(tile1_seq.split(' ')[1]))
              bad_conflicts += str('{0:.2f}'.format(eval_colocated_end_pair(tile2_seq.split(' ')[0], tile1_seq.split(' ')[2], temp_in_C)))+' = eval_colocated_end_pair('+str(tile2_seq.split(' ')[0])+','+str(tile1_seq.split(' ')[2])+','+str(temp_in_C)+')\n'
      
  if not threaded:
    with open(directory+'bad_arbitrary_lattice_conflicts - '+lattice_spacer+'.txt', 'w') as f:
        f.write('Showing erroneous attachments with energy $\leq ' +str(thres)                     +'$ kcal/mol for single-match and single-mismatch binding events '                     +'to arbitrary (correct or incorrect) lattices. These errors '                     +'may or may not respect proofreading block tile-position.')
        f.write('\\begin{footnotesize}\\begin{verbatim}\n')
        f.write(bad_conflicts)
        f.write('\\end{verbatim}\\end{footnotesize}\n')
                  
  histogram([correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies],
            labels=['correct growth ({} energies)'.format(len(correct_growth_energies)),
                    'one mismatch below ({} energies)'.format(len(single_mismatch_below_energies)),
                    'one mismatch above ({} energies)'.format(len(single_mismatch_above_energies))
                   ],
            xaxis='energy (-kcal/mol)',yaxis='',
            label_size=12,legend_font_size=9,
            title='growth on right (detailed binding, incorrect and correct lattices, lattice spacer '+lattice_spacer*(bool(lattice_spacer))+'is empty string'*(not(bool(lattice_spacer)))+')', #; for {} tiles)'.format(len(names_seqs)),
            plot_dir=directory,xmin=-18,num_bins=num_bins) # xmin=4,xmax=18
  
  sys.stdout.write('\n')
  return correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies


def lattice_binding_energies_simplified(f,directory, temp_in_C=53, lattice_spacer='TTTTT',num_bins=-1):
  '''Nupack binding() energy of tile binding of right of nanotube with no mistmatches
  or one mismatch'''    
  names, _, _, seqs, seqs_bt  = read_strand_data_from_idt_order(f)
  domains, domains_incl_biotin_labels = read_domains_from_idt_order(f)
  if not os.path.exists(directory): os.makedirs(directory)
  print('lattice_binding_energies_simplified(), lattice binding analysis from file:\n'+f+' plots will be placed in: ' + directory)
  names_seqs = list(zip(names,seqs))
   
  threaded = _threaded
  num_tiles_processed = 0
  correct_growth_energies = []
  single_mismatch_below_energies = []
  single_mismatch_above_energies = []

  for tile1_name, tile1_seq in names_seqs:
    tile1_domains = tile1_seq.split(' ')
    sys.stdout.write('.')
    if num_tiles_processed%10==0 and num_tiles_processed>0: sys.stdout.write(str(num_tiles_processed))
    #if num_tiles_processed%10==0: sys.stdout.write('.')
    sys.stdout.flush()
    num_tiles_processed += 1

    #correct_growth_energies.append(binding(tile1_domains[1],wc(tile1_domains[1]),temp_in_C)+binding(tile1_domains[2],wc(tile1_domains[2]),temp_in_C))
    correct_growth_energies.append(binding(tile1_domains[1]+tile1_domains[2], 
                                           wc(tile1_domains[2])+lattice_spacer+wc(tile1_domains[1]),
                                           temp_in_C))
    
    # mismatch below
    results1 = [global_thread_pool.apply_async(binding, args=(tile1_domains[1]+tile1_domains[2],
                wc(tile1_domains[2])+lattice_spacer+tile2_seq.split(' ')[3],
                temp_in_C)) for _,tile2_seq in names_seqs if wc(tile1_seq.split(' ')[1])!=tile2_seq.split(' ')[3]]
    
    single_mismatch_below_energies = single_mismatch_below_energies + [result.get() for result in results1] 

    # mismatch above
    results2 = [global_thread_pool.apply_async(binding, args=(tile1_domains[1]+tile1_domains[2],
                tile2_seq.split(' ')[0]+lattice_spacer+wc(tile1_domains[1]),
                temp_in_C)) for _,tile2_seq in names_seqs if wc(tile1_seq.split(' ')[2])!=tile2_seq.split(' ')[0] ]
    
    single_mismatch_above_energies = single_mismatch_above_energies + [result.get() for result in results2] 

  histogram([correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies],
            labels=['correct growth ({} energies)'.format(len(correct_growth_energies)),
                    'one mismatch below ({} energies)'.format(len(single_mismatch_below_energies)),
                    'one mismatch above ({} energies)'.format(len(single_mismatch_above_energies))
                   ],
            xaxis='energy (-kcal/mol)',yaxis='',
            title='growth on right (simplified binding(), lattice spacer is '+lattice_spacer*(bool(lattice_spacer))+'empty string'*(not(bool(lattice_spacer)))+')', 
            label_size=12,legend_font_size=9,
            plot_dir=directory,xmin=-18,num_bins=num_bins) # ,xmin=4,xmax=18

  return correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies
 


def tile_helix(tile_name):
  '''return helix number for tile with name tile_name'''
  n = tile_num(tile_name)
  d = tile_dir(tile_name) 
  #if n%2==0:
  if d == 'nw' or d == 'se' : return 2*n
  elif d == 'sw': return 2*n-1
  elif d == 'ne': return (2*n+1) % _num_helices
  else: sys.exit('Exiting program: Error in tile name. Was expecting a tile name with two semi-colons followed by one of nw,ne,se,sw. E.g. U6;00->00;ne.')

def tile_dir(name):
  return name[10:12]
def tile_num(name):
  return int(name[1:2])


def analyse_algo_conflicts(f,directory,temp_in_C=53,lattice_spacer='TTTTT',threaded=True,num_bins=-1):
  '''binding() energy of algorithmic conflicting tile (single-mismatches that can arise 
  as the first error from a locally-correct lattice)'''

  names, _, _, seqs, seqs_bt  = read_strand_data_from_idt_order(f)
  domains, domains_incl_biotin_labels = read_domains_from_idt_order(f)
  
  if not os.path.exists(directory): os.makedirs(directory)
  print('  the next few plots will be placed in: ' + directory)
  names_seqs = list(zip(names,seqs))
  print('  analysing single mismatches that are the first error in a \
locally correct lattice (aka first algorithmic conflicts) using \
the function analyse_algo_conflicts(), from file: ' + f +\
' (iterating through '+str(len(names_seqs))+' strands):\n')  
  
  num_tiles_processed = 0
  correct_growth_energies = []
  single_mismatch_below_energies = []
  single_mismatch_above_energies = []

  thres = _bad_first_algo_conflict_threshold
  lower_thres_correct_binding = _lower_thres_correct_binding
  bad_conflicts = ''
  weak_correct_bindings = ''
       
  for tile1_name, tile1_seq in names_seqs:
    tile1_domains = tile1_seq.split(' ')
    if threaded:    
      if num_tiles_processed%10==0 and num_tiles_processed>0: sys.stdout.write(str(num_tiles_processed))
      sys.stdout.write('.')
      sys.stdout.flush()
      num_tiles_processed += 1

      correct_growth_energies.append(lattice_binding_spacer(tile1_seq,
                                                     wc(tile1_seq.split(' ')[2]),
                                                     lattice_spacer,
                                                     wc(tile1_seq.split(' ')[1]),
                                                     temp_in_C))
      
      mismatch_below_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
                                                           wc(tile1_seq.split(' ')[2]),
                                                           lattice_spacer,
                                                           wc(tile2_seq.split(' ')[1]),
                                                           temp_in_C)) \
                            for _,tile2_seq in names_seqs \
                                if tile1_seq.split(' ')[2]==tile2_seq.split(' ')[2] and \
                                   tile1_seq.split(' ')[1]!=tile2_seq.split(' ')[1]]
      single_mismatch_below_energies = single_mismatch_below_energies + [result.get() for result in mismatch_below_res] 

      mismatch_above_res = [global_thread_pool.apply_async(lattice_binding_spacer, args=(tile1_seq,
                                                           wc(tile2_seq.split(' ')[2]),
                                                           lattice_spacer,
                                                           wc(tile1_seq.split(' ')[1]),
                                                           temp_in_C)) \
                            for _,tile2_seq in names_seqs \
                            if tile1_seq.split(' ')[1]==tile2_seq.split(' ')[1] and \
                               tile1_seq.split(' ')[2]!=tile2_seq.split(' ')[2]]
      single_mismatch_above_energies = single_mismatch_above_energies + [result.get() for result in mismatch_above_res] 

    else: # not threaded
      e = lattice_binding_spacer(tile1_seq,
                                  wc(tile1_seq.split(' ')[2]),
                                  lattice_spacer,
                                  wc(tile1_seq.split(' ')[1]),
                                  temp_in_C)
      correct_growth_energies.append(e)
      if e > lower_thres_correct_binding:
        weak_correct_bindings += print_lattice_binding(e, tile1_seq, tile1_name)
      
      for tile2_name,tile2_seq in names_seqs:
        # mismatch below
        if tile1_seq.split(' ')[2]==tile2_seq.split(' ')[2] and                    tile1_seq.split(' ')[1]!=tile2_seq.split(' ')[1]:
          e = lattice_binding_spacer(tile1_seq,
                              wc(tile1_seq.split(' ')[2]),
                              lattice_spacer,
                              wc(tile2_seq.split(' ')[1]), # there was previously an eror here
                              temp_in_C)
          single_mismatch_below_energies.append(e)
          if e < thres: 
              #print_lattice_mismatch_below(energy, top_output, tile_seq, tile_name, bottom_output, bottom_name):
            bad_conflicts += print_lattice_mismatch_below(e, wc(tile1_seq.split(' ')[2]), tile1_seq, tile1_name, wc(tile2_seq.split(' ')[1]), 'mismatch: dom 3 of tile WC to dom 1 of '+tile2_name)
            bad_conflicts += str('{0:.2f}'.format(eval_colocated_end_pair(tile1_seq.split(' ')[1],wc(tile2_seq.split(' ')[1]),temp_in_C))) +' = eval_colocated_end_pair('+str( tile1_seq.split(' ')[1])+','+str(wc(tile2_seq.split(' ')[1]))+','+str(temp_in_C)+')\n'
        # mismatch above
        if tile1_seq.split(' ')[1]==tile2_seq.split(' ')[1] and                    tile1_seq.split(' ')[2]!=tile2_seq.split(' ')[2]:
          e = lattice_binding_spacer(tile1_seq,
                                wc(tile2_seq.split(' ')[2]),
                                lattice_spacer,
                                wc(tile1_seq.split(' ')[1]),
                                temp_in_C)
          single_mismatch_above_energies.append(e)
          if e < thres:
            bad_conflicts += print_lattice_mismatch_above(e, wc(tile2_seq.split(' ')[2]), 'mismatch: dom 0 of tile WC to dom 2 of '+tile2_name, tile1_seq, tile1_name, wc(tile1_seq.split(' ')[1]))
            bad_conflicts += str('{0:.2f}'.format(eval_colocated_end_pair(wc(tile2_seq.split(' ')[2]), tile1_seq.split(' ')[2], temp_in_C)))+' = eval_colocated_end_pair('+str(wc(tile2_seq.split(' ')[2]))+','+str(tile1_seq.split(' ')[2])+','+str(temp_in_C)+')\n'

  if not threaded:
    with open(directory+'bad_first_algo_conflicts - '+lattice_spacer+'.txt', 'w') as f:
        f.write('Showing erroneous attachments with energy $\leq ' +str(thres)                     +'$ kcal/mol for single mismatches (algorithmic conflicts) that occur as '                     +'the first error in a proofreading block.')
        f.write('\\begin{footnotesize}\\begin{verbatim}\n')
        f.write(bad_conflicts)
        f.write('\\end{verbatim}\\end{footnotesize}\n')
        #f.write('Done writing ' +str(len(bad_conflicts))+ ' characters')
        #f.close()
    with open(directory+'weak_correct_bindings - '+lattice_spacer+'.txt', 'w') as f:
        f.write('Showing correct lattice attachments with energy $\geq '+str(lower_thres_correct_binding)+'$ kcal/mol.')
        f.write('\\begin{footnotesize}\\begin{verbatim}\n')
        f.write(weak_correct_bindings)
        f.write('\\end{verbatim}\\end{footnotesize}\n')
  

  histogram([correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies],
            labels=['correct growth ({} energies)'.format(len(correct_growth_energies)),
                    'one mismatch below ({} energies)'.format(len(single_mismatch_below_energies)),
                    'one mismatch above ({} energies)'.format(len(single_mismatch_above_energies))
                   ],
            xaxis='energy (kcal/mol)',yaxis='count',
            title='correct and algorithmic error attachments to a valid lattice'+' (lattice spacer is empty)'*(not(bool(lattice_spacer))), #; for {} tiles)'.format(len(names_seqs)),
            filename='algorithmic error (first), dG before and after, mismatches above and below, lattice spacer is '+lattice_spacer*(bool(lattice_spacer))+'empty string'*(not(bool(lattice_spacer))), #; for {} tiles)'.format(len(names_seqs)),            
            normed=0,label_size=12,legend_font_size=9,plot_dir=directory,
            xmin=min(-17.1,min(correct_growth_energies+single_mismatch_below_energies+single_mismatch_above_energies)),
            xmax=0,
            num_bins=num_bins) # was ymax=282
  histogram([correct_growth_energies,single_mismatch_below_energies+single_mismatch_above_energies],
            labels=['{} correct attachments'.format(len(correct_growth_energies)),
                    '{} algorithmic errors'.format(len(single_mismatch_below_energies+single_mismatch_above_energies))
                   ],
            xaxis='energy (kcal/mol)',yaxis='count',
            title=   'correct and algorithmic error attachments to a valid lattice'+' (lattice spacer is empty)'*(not(bool(lattice_spacer))), #; for {} tiles)'.format(len(names_seqs)),
            filename='algorithmic error (first), dG before and after, lattice spacer is '+lattice_spacer*(bool(lattice_spacer))+'empty string'*(not(bool(lattice_spacer)))+')', #; for {} tiles)'.format(len(names_seqs)),
            normed=0,label_size=12,legend_font_size=9,plot_dir=directory,
            xmin=min(-17.1,min(correct_growth_energies+single_mismatch_below_energies+single_mismatch_above_energies)),
            xmax=0,
            num_bins=num_bins) # was ymax=282

  sys.stdout.write('\n')
  return correct_growth_energies,single_mismatch_below_energies,single_mismatch_above_energies



def print_lattice_binding(energy, tile_seq, tile_name):
  s = '\n'
  s += str('{0:.2f}'.format(energy)) + ' binding strength\n'                        
  s += '<-' + wc(tile_seq.split(' ')[2])[::-1]+'\n'
  s += ' /' +tile_seq.split(' ')[2]+'-'+tile_seq.split(' ')[3]+'->  '+tile_name+'\n'
  s += ' \\' +tile_seq.split(' ')[1][::-1]+'-'+tile_seq.split(' ')[0][::-1]+'\n'
  s += ' -' + wc(tile_seq.split(' ')[1]) +'->\n'
  return s

def print_lattice_mismatch_below(energy, top_output, tile_seq, tile_name, bottom_output, bottom_name):
    s = '\n'
    s += str('{0:.2f}'.format(energy)) + ' strength mismatch\n'                        
    s += '<-' + top_output[::-1]+'\n'
    s += ' /' +tile_seq.split(' ')[2]+'-'+tile_seq.split(' ')[3]+'->  '+tile_name+'\n'
    s += ' \\' +tile_seq.split(' ')[1][::-1]+'-'+tile_seq.split(' ')[0][::-1]+'\n'
    s += ' -' + bottom_output +'->                '+bottom_name +'\n'
    #print s ; sys.stdout.flush()
    return s

def print_lattice_mismatch_above(energy, top_output, top_name, tile_seq, tile_name, bottom_output):
    s = '\n'
    s += str('{0:.2f}'.format(energy)) + ' strength mismatch\n'                        
    s += '<-' + top_output[::-1]+'                '+top_name + '\n'
    s += ' /' +tile_seq.split(' ')[2]+'-'+tile_seq.split(' ')[3]+'->  '+tile_name+'\n'
    s += ' \\' +tile_seq.split(' ')[1][::-1]+'-'+tile_seq.split(' ')[0][::-1]+'\n'
    s += ' -' + bottom_output +'->\n'
    #print s; sys.stdout.flush()
    return s

def analyse_domains(names,seqs,seqs_bt,domains,domains_incl_biotin_labels,
                    seq_ws,seq_ws_bt,temp_in_C=53.0,num_domain_pairs=False,
                    directory='plots_tmp/',non_standard_tiles=0,pickle_data=0,num_bins=-1):
  '''Analyse domain sequences'''
  
  print('\ndomain analysis: plots will be placed in: ' + directory)
  #print 'Analysing sequences from file:\n' + filename + '\n' +str(len(domains))+' unique domains, ' +str(len(domains_incl_biotin_labels))+' domains (leaving in biotin labels), ' + str(len([d for d in domains if len(d)==10])) + ' of length 10, '+str(len([d for d in domains if len(d)==11])) + ' of length 11' ; sys.stdout.flush()

  print('  analysing tile input and output sec structure'); sys.stdout.flush()

  if non_standard_tiles:
    print(seq_ws)
    seq_ws = [s.split()[0]+' '+s.split()[2]+' '+s.split()[3]+' '+s.split()[5] for s in seq_ws if len(s.split())==6]
    print(seq_ws)

  input_energies = analyse_input_pairs(seq_ws, temp_in_C)
  output_energies = analyse_output_pairs(seq_ws, temp_in_C)   
  lattice_input_energies_correct_binding = analyse_lattice_input_energies_correct_binding(seq_ws, temp_in_C)
  histogram([input_energies,lattice_input_energies_correct_binding,output_energies],
            labels=['input','lattice (correct)','output'],xaxis='energy (kcal/mol)',yaxis='',
            title='tile input, lattice and output secondrary structure (ss)',plot_dir=directory,num_bins=num_bins)
  scatter_plot_energies(input_energies, lattice_input_energies_correct_binding, xaxis='input sec struct', yaxis='lattice sec struct (correct binding)', title_info='',plot_dir=directory) #,xthreshold=9.3,ythreshold=9.3
  
  # domain secondary structure
  sys.stdout.write('  analysing domain secondary structure: '); sys.stdout.flush()
  domain_energies = [pfunc(d, temp_in_C) for d in domains] 
  sys.stdout.write(str(len(domain_energies))+ ' domains\n')

  strong_doms = [d.replace('/iBiodT/','T') for d in domains_incl_biotin_labels if '/iBiodT/' in d]
  # in the following line I should replace d.replace('/iBiodT/','T') with d:
  weak_doms = [d for d in domains_incl_biotin_labels if '/iBiodT/' not in d]    
  strong_doms_energies = [pfunc(d, temp_in_C) for d in strong_doms] 
  weak_doms_energies = [pfunc(d, temp_in_C) for d in weak_doms] 
  if len(strong_doms_energies)>0 and len(weak_doms_energies):
    histogram([weak_doms_energies,strong_doms_energies],labels=['domains','domains with biotin'],
              xaxis='energy (kcal/mol)',yaxis='',title='domain secondary structure energy histogram (via NUPACK pfunc)',
              plot_dir=directory,num_bins=num_bins)
  else:
    #weak_doms_energies.extend(strong_doms_energies)
    histogram([domain_energies],labels=['domains'],xaxis='energy (kcal/mol)',yaxis='',
      title='domain secondary structure (pfunc)',plot_dir=directory,num_bins=num_bins)
  
  if num_domain_pairs != -1:
    sys.stdout.write('  analysing domain pairs: ');sys.stdout.flush()
    domains_pairs = get_seq_pairs(domains)
    print(str(len(domains_pairs)) +  ' domain pairs to process'); sys.stdout.flush()

    if num_domain_pairs>0: 
      domains_pairs = domains_pairs[:num_domain_pairs]
      print('Warning: I am merely going to analyse a strict subset ('+str(len(domains_pairs))+') of all domain pairs')
  
  #if num_domain_pairs: domains_pairs = domains_pairs[:num_domain_pairs]
  #sys.stdout.write('  analysing domain pair secondary structure: '); sys.stdout.flush()
  #print(str(len(domains_pairs)) +  ' strand pairs to process');sys.stdout.flush()
  
    non_wc_domains_pairs = [p for p in domains_pairs if p[0]!=wc(p[1])]
    wc_domains_pairs = [p for p in domains_pairs if p[0]==wc(p[1])]
    RNAduplex_non_wc_domains_pairs_energies = tile_pair_RNAduplex_energies(non_wc_domains_pairs,temp_in_C)
    # unpaired strands get RNAduplex energy of ~ -100, so set those to -5:
    if MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE:
      RNAduplex_non_wc_domains_pairs_energies = [e if e > -100 else -5 for e in RNAduplex_non_wc_domains_pairs_energies ]
    else:
      RNAduplex_non_wc_domains_pairs_energies = [e if e < 100 else 5 for e in RNAduplex_non_wc_domains_pairs_energies ]

    threshold = 3.0
  
    print(str(len( [e for e in RNAduplex_non_wc_domains_pairs_energies if e > threshold ]))+' non-wc domain pairs with energy > '+str(threshold))
  
    RNAduplex_wc_domains_pairs_energies = tile_pair_RNAduplex_energies(wc_domains_pairs,temp_in_C)

    histogram([RNAduplex_wc_domains_pairs_energies,RNAduplex_non_wc_domains_pairs_energies],
            labels=['wc domains','non-wc domains'],xaxis='energy (kcal/mol)',yaxis='',
            title='domain pair (RNAduplex)',plot_dir=directory,num_bins=num_bins)
  
    non_wc_binding_pair_energies = tile_pair_NUPACK_binding_energies(non_wc_domains_pairs,temp_in_C)
    scatter_plot_energies(non_wc_binding_pair_energies, RNAduplex_non_wc_domains_pairs_energies, xaxis='nupack domain pair binding', yaxis='RNAduplex domain pair (1999)', title_info='non-wc',plot_dir=directory) #,xthreshold=9.3,ythreshold=9.3

    wc_binding_pair_energies = tile_pair_NUPACK_binding_energies(wc_domains_pairs,temp_in_C)
    scatter_plot_energies(wc_binding_pair_energies, RNAduplex_wc_domains_pairs_energies, xaxis='nupack domain pair binding', yaxis='RNAduplex domain pair (1999)', title_info='wc',plot_dir=directory) #,xthreshold=9.3,ythreshold=9.3

    histogram([wc_binding_pair_energies,non_wc_binding_pair_energies],
            labels=[str(len(wc_binding_pair_energies))+ ' WC domain pairs',str(len(non_wc_binding_pair_energies))+' non-WC domain pairs'],
            xaxis='energy (kcal/mol)',yaxis='',
            normed=1,label_size=12,legend_font_size=9,
            title='domain pair (nupack binding())',plot_dir=directory,num_bins=num_bins)





#############################################################################
# Reading from idt-formatted files (with extension .idt) 
#############################################################################

def get_lines_from_file(filename):
    '''Reads uncommented and non-empty lines from an IDT-formatted file'''
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines = [line[:line.find('#')].strip() for line in lines]
    lines = [line for line in lines if len(line) > 0]
    return lines

def is_DNA_sequence(s):
    return set(s.replace('/iBiodT/','').replace(' ','')).issubset('ACTG')

def get_DNA_sequences_from_idt_file(f):
    return get_testtube_DNA_sequences_from_idt_file(f) +  get_plate_DNA_sequences_from_idt_file(f)

def get_testtube_lines_from_idt_file(f):
    lines = get_lines_from_file(f)
    return [l for l in lines if len(l.split(','))==4 and is_DNA_sequence(l.split(',')[1]) and 'guard' not in l]
        
def get_plate_lines_from_idt_file(f):
    '''Ignores any line containing the word guard'''
    lines = get_lines_from_file(f)
    return [l for l in lines if len(l.split(','))==4 and is_DNA_sequence(l.split(',')[2]) and 'guard' not in l]
        
def read_strand_data_from_idt_order(f):
    tiles = get_testtube_lines_from_idt_file(f) 
    tile_names = [l.split(',')[0].replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for l in tiles]
    tile_seqs_without_biotin_labels = [tile.split(',')[1].replace('/iBiodT/','T').replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]
    tile_seqs_with_biotin_labels = [tile.split(',')[1].replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles] 
    tile_seqs_with_whitespace_between_doms = [tile.split(',')[1].replace('/iBiodT/','T').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]
    tile_seqs_biotins_whitespace_between_doms = [tile.split(',')[1].replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]    
    tiles = get_plate_lines_from_idt_file(f) 
    tile_names += [l.split(',')[1].replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for l in tiles]
    tile_seqs_without_biotin_labels += [tile.split(',')[2].replace('/iBiodT/','T').replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]
    tile_seqs_with_biotin_labels += [tile.split(',')[2].replace(' ','').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles] 
    tile_seqs_with_whitespace_between_doms += [tile.split(',')[2].replace('/iBiodT/','T').replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]
    tile_seqs_biotins_whitespace_between_doms += [tile.split(',')[2].replace('\n','').replace('\r','').lstrip().strip() for tile in tiles]    
    
    if len(tile_seqs_without_biotin_labels) == 0: 
      exit('There are 0 DNA strands. Did you give the correct input? Strands should be specified in idt file format. Run with --help to get an example.')
    if len(tile_names) != len(tile_seqs_without_biotin_labels): 
      exit('Error: there are '+str(len(tile_names))+ ' strand names and '+str(len(tile_seqs_without_biotin_labels))+'strands!')
    return tile_names, tile_seqs_without_biotin_labels, tile_seqs_with_biotin_labels, tile_seqs_with_whitespace_between_doms, tile_seqs_biotins_whitespace_between_doms

    
def read_domains_from_idt_order(filename, non_standard_tiles=0):
    _, _, _, tile_seqs_with_whitespace_between_doms, tile_seqs_biotins_whitespace_between_doms = read_strand_data_from_idt_order(filename)
    doms = []
    doms_with_biotin = []
    for l in tile_seqs_with_whitespace_between_doms:
        doms.extend(l.replace('\n','').replace('\r','').split(' ') )
        
    for l in tile_seqs_biotins_whitespace_between_doms:    
        doms_with_biotin.extend(l.replace('\n','').replace('\r','').split(' ') )
        
    return list(set(doms)), list(set(doms_with_biotin))


def flatten(l):
  return [item for sublist in l for item in sublist]


#############################################################################
# Plots 
#############################################################################

def histogram(data, labels='', xaxis='', yaxis='', title='', filename='', plot_dir='', log_scale=False, 
              normed=1, xmin=0, xmax=0, ymin=0, ymax=0, label_size=12, legend_font_size=12, show_mean_min_max=1, 
              legend_location='upper left',num_bins=-1):
  '''PLot a histogram. By default data is normalized (area under curve 
  sums to 1 for each dataset (curve colour))'''
  min_data = 10000; max_data = -10000

  # Fix up function param types/formatting
  if filename=='': filename=title
  if type(data)!=list: data = [data]
  if labels=='':labels=['' for _ in data]
  if title!='': filename = title
  else: filename = '_'.join(l for l in labels)
      
  pp = PdfPages(plot_dir + filename + '.pdf')
  colours = itertools.cycle(['r','b','g','y'])
  
  for i,d in enumerate(data):
    if d!= []:     #assert len(d)!=0
      min_data = min(min_data,min(d)) 
      max_data = max(max_data,max(d)) 

  if num_bins==-1: num_bins = min(50, math.ceil(math.sqrt(len(flatten(data))) +1) )

  step_size = (max_data-min_data)/num_bins
  bins = [min_data+(i*step_size) for i in range(num_bins)  ]

  for i,d in enumerate(data):
    if d != []:
      # ax = (y,x,_) is a 
      ax = plt.hist(d, bins, normed=normed, 
                  range=[2.0,18.0],
                  histtype='stepfilled', 
                  facecolor=next(colours), alpha=0.5, edgecolor='black', 
                  label=labels[i] + show_mean_min_max*((labels[i]!='')*', '+'mean='+str('{0:.2f}'.format(sum(d)/float(len(d))))+', min='+str('{0:.2f}'.format(min(d)))+', max='+str('{0:.2f}'.format(max(d)))))
  
  if xmin!=0 or xmax!=0: plt.xlim(xmin, xmax)
  if ymin!=0 or ymax!=0: 
    plt.ylim(ymin, ymax)
  else:
    plt.ylim(ax[0].min(),  ax[0].max()*1.4  )

  plt.legend(loc=legend_location, fontsize=legend_font_size) #'smaller'
  plt.title(title, fontsize = label_size) #'smaller'
  plt.ylabel(yaxis,fontsize = label_size)
  plt.xlabel(xaxis,fontsize = label_size)
  # Set the tick labels font
  
  plt.xticks(fontsize = label_size)
  plt.yticks(fontsize = label_size)
  plt.tick_params(axis='y',which='minor',bottom='off')
  
  plt.minorticks_on()
  if log_scale: plt.yscale('log', nonposy='clip')
  pp.savefig()
  plt.clf()
  pp.close()

  
def scatter_plot_energies(xaxis_energies,yaxis_energies,xaxis='',yaxis='',xthreshold=False,ythreshold=False,title_info='',plot_dir='',label_size=12):
    '''Generate a scatter plot comparing two equal-length lists of floats'''
    
    if len(xaxis_energies) != len(yaxis_energies):
        print('Unequal length energy lists in scatter_plot_energies(...): ' + str(len(xaxis_energies)) +' != '+str(len(yaxis_energies)))
        return None
    
    #title = xaxis + ' V ' + yaxis + ' - ' + str(len(xaxis_energies)) + ' tile pairs' + (' - '+title_info  if title_info != '' else '')  
    title = xaxis + ' V ' + yaxis + (' - '+title_info  if title_info != '' else '')  
    filename = title
    # if threshold != False: title += ' thres='+str(threshold)
    timestamp = datetime.fromtimestamp(time.time()).strftime('%m-%d-%H-%M-%S')
    pp = PdfPages(plot_dir + filename + '.pdf')
        #'tilesets/plots/' + filename + ' - ' + timestamp + '.pdf')
    

    marker='.'
    if len(xaxis_energies) < 101: marker='o'
    ps = plt.scatter(xaxis_energies,yaxis_energies, 
                    # c = my_colours, cmap=Cmap,   # , cmap='gray') # ,s=0.5
                     alpha=0.5, marker=marker,edgecolors='none') 
    plt.grid(True)
    
    plt.annotate('min '+xaxis+' = '+"%.1f" % min(xaxis_energies), 
                 xy=(.45, .07), xycoords='axes fraction')# xy=(7,min(yaxis_energies)+0.5))
    plt.annotate('min '+yaxis+' = '+"%.1f" % min(yaxis_energies), 
                 xy=(.45, .12), xycoords='axes fraction') # xy=(7,min(yaxis_energies)+1))
    if xthreshold != False:
        # plt.plot([2.1, 14], [threshold, threshold], 'k-')
        plt.annotate(str(count_pairs_above_threshold(xaxis_energies,xthreshold))+' of '+str(len(xaxis_energies))+' with '+xaxis+' > '+str(xthreshold), 
                     xy=(.45, .15), xycoords='axes fraction') #xy=(7,min(yaxis_energies)+1.5) ) 
    if ythreshold != False:
        # plt.plot([2.1, 14], [threshold, threshold], 'k-')
        #above = count_pairs_above_threshold(yaxis_energies, threshold)   
        plt.annotate(str(count_pairs_above_threshold(yaxis_energies,ythreshold))+' of '+str(len(yaxis_energies)) +' with '+yaxis+' > '+str(ythreshold), 
                     xy=(.45, .2), xycoords='axes fraction')  # xy=(7,min(yaxis_energies)+2) ) 

    plt.title(title,fontsize=label_size)
    plt.ylabel(yaxis,fontsize=label_size)
    plt.xlabel(xaxis,fontsize=label_size)
    plt.tight_layout()
    plt.minorticks_on()
    pp.savefig()
    plt.clf()
    pp.close()

def count_pairs_above_threshold(eng, t):
    count = 0
    for e in eng:
        if e > float(t)  : count +=1
    return count


def plot_energies2(list_of_pairs_of_lists, labels=[''], xaxis='', yaxis='', title_info='',plot_dir='',label_size=12):
    for pair in list_of_pairs_of_lists:
        if len(pair[0]) != len(pair[1]):
            print('Error: Unequal length energy lists in function plot_energies2(..): ' + str(len(pair[0])) +' != '+str(len(pair[1])))
            return None
    
    # title = xaxis + ' V ' + yaxis + ' - ' + str(len(list_of_pairs_of_lists[0][0])) + ' tile pairs' + (' - '+title_info  if title_info != '' else '')  
    title = xaxis + ' V ' + yaxis + (' - '+title_info  if title_info != '' else '')  
    filename = title
    # if threshold != False: title += ' thres='+str(threshold)
    timestamp = datetime.fromtimestamp(time.time()).strftime('%m-%d-%H-%M-%S')
    pp = PdfPages(plot_dir + filename + '.pdf')
        #'tilesets/plots/' + filename + ' - ' + timestamp + '.pdf')
    colours = itertools.cycle(['r','b','g','y'])
    fig, ax = plt.subplots()
    
    marker='.'
    if len(list_of_pairs_of_lists) < 101: marker='o'
    for i,pair in enumerate(list_of_pairs_of_lists):
        plt.scatter(pair[0],pair[1], 
                     alpha=0.5, marker=marker,edgecolors='none',color=next(colours), label=str(labels[i])+' ('+str(len(pair[0]))+' values)' ) 
    
    #start, end = ax.get_xlim()
    #ax.xaxis.set_ticks(np.arange(start, end, 2.0))
    plt.grid(True)
    
    plt.legend(loc='lower right')
    plt.title(title,fontsize=label_size)
    plt.ylabel(yaxis,fontsize=label_size)
    plt.xlabel(xaxis,fontsize=label_size)
    plt.tight_layout()
    plt.minorticks_on()
    pp.savefig()
    plt.clf()
    pp.close()

def convert_scatter_plots_to_jpg(f, directory):
  # This function is causing problems: perhaps due to version mismatch between imagemagik and wand
  directory = directory 
  print('Converting some plots from pdf to jpg. Working in:\n' + directory) ; sys.stdout.flush()
  #directory = _root_dir_analysis + f[1]+'/'     
  for pdf_file in filelist:
    if os.path.isfile(directory+pdf_file):
      print('Converting pdf to jpg file:\n    '+directory+pdf_file +'\n--> '+directory+pdf_file[:-4]+".jpg\n")
      with Image(filename=directory+pdf_file) as img:
        img.save(filename=directory+pdf_file[:-4]+".jpg")
  print('done converting pdf images to jpg images')



#############################################################################
# Run analysis on analyse tile pairs 
#############################################################################


def get_time():
    return time.time()

def split(l, num):
    """Yield successive chunks from l, each of same size, to make num total.    
    http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python"""
    assert 0 < num <= len(l)
    n = len(l) // num
    if len(l) % num != 0:
        n += 1
    for i in range(0, len(l), n):
        yield l[i:i+n]

def get_seq_pairs(seqs):
    seqs_pairs = []
    for i, str_i in enumerate(seqs):
        for j, str_j in enumerate(seqs[i+1:]):
            seqs_pairs.append( [ str_i, str_j ] )
    return seqs_pairs

def tile_pair_RNAduplex_energies(sp,temp_in_C,NA_parameter_set='dna_mathews1999.par'):
    sys.stdout.write('  starting RNAduplex_multiple tile-pair job  (using parameter set '+str(NA_parameter_set)+') on ' + str(len(sp))+ ' strand pairs at time: '); start = get_time(); sys.stdout.write( str(datetime.now().time()) )
    if 0 < global_thread_pool._processes < len(sp): #_threaded:
        list_of_list_of_seqpair = list(split(sp, global_thread_pool._processes))
        sys.stdout.write(' (processing %d lists of sequence pairs, each of size ~%d);' % (len(list_of_list_of_seqpair), len(list_of_list_of_seqpair[0])))
        results = [global_thread_pool.apply_async(RNAduplex_multiple, args=(s, temp_in_C, NA_parameter_set)) for s in list_of_list_of_seqpair]
        list_of_list_of_energies = [result.get() for result in results]
        e = [energy for list_of_energies in list_of_list_of_energies for energy in list_of_energies]
    else:
        e = RNAduplex_multiple(sp, temp_in_C, NA_parameter_set)
    end = get_time(); print(' duration: %.2f seconds, for %d pairs ' % ((end-start), len(e))); sys.stdout.flush()
    return e

#print '  processing %.2f lists of sequence pairs, each of size ~%d' % (len(list_of_list_of_seqpair), len(list_of_list_of_seqpair[0]))
def tile_pair_NUPACK_binding_energies(sp,temp_in_C):
    sys.stdout.write('  starting NUPACK strand-pair binding job on ' + str(len(sp))+ ' strand pairs at time: '); start = get_time(); sys.stdout.write(str(datetime.now().time())+';'); sys.stdout.flush(); 
    if _threaded:
        results = [global_thread_pool.apply_async(binding, args=(p[0], p[1], temp_in_C)) for p in sp]
        e = [result.get() for result in results]
    else:
        e = [ binding(p[0],p[1],temp_in_C )  for p in sp]  
    end = get_time(); print(' duration: %.2f seconds, for %d pairs ' % ((end-start), len(e))); sys.stdout.flush()
    return e

def tile_pair_NUPACK_mfe_energies(sp,temp_in_C):
    sys.stdout.write('  starting NUPACK tile-pair mfe job on ' + str(len(sp))+ ' tile pairs at time: '); start = get_time(); sys.stdout.write(str(start)+';'); sys.stdout.flush(); 
    if _threaded and False: #  threading with vanilla mfe() causes problems due to files overwriting each other  getting en arror so turning threaded off: mfe uses file I/O which does not work with 'threaded'  
        results = [global_thread_pool.apply_async(mfe_binding_async, args=(p[0], p[1], temp_in_C)) for p in sp]
        e = [result.get() for result in results]
    else:
        e = [ mfe_binding(p[0],p[1],temp_in_C)  for p in sp]  
    end = get_time(); print(' duration: %.2f seconds, for %d pairs ' % ((end-start), len(e))); sys.stdout.flush()
    return e



def analyse_input_pairs(seq, temp_in_C):
    seqs_4_domains = [s for s in seq if len(s.split(' '))==4]
    print(("  analysing input pair energy for "+str(len(seqs_4_domains))+" sequences with 4 domains; ignoring "+str(len(seq)-len(seqs_4_domains))+" sequences that do not have 4 domains."))
    results = [global_thread_pool.apply_async(pfunc, args=(s.split(' ')[1]+s.split(' ')[2],temp_in_C)) for s in seqs_4_domains]
    energies = [result.get() for result in results]
    return energies

def analyse_output_pairs(seq_ws, temp_in_C):
    seqs_ws_4_domains = [s for s in seq_ws if len(s.split(' '))==4]
    print(("  analysing output pair energy for "+str(len(seqs_ws_4_domains))+" sequences with 4 domains; ignoring "+str(len(seq_ws)-len(seqs_ws_4_domains))+" sequences that do not have 4 domains."))
    results = [global_thread_pool.apply_async(pfunc, args=(s.split(' ')[0]+s.split(' ')[3],temp_in_C)) for s in seqs_ws_4_domains]
    energies = [result.get() for result in results]    
    return energies

def analyse_lattice_input_energies_correct_binding(seq_ws, temp_in_C):
    seqs_ws_4_domains = [s for s in seq_ws if len(s.split(' '))==4]
    print(("  analysing lattice_input_energies_correct_binding pair energy for "+str(len(seqs_ws_4_domains))+" sequences with 4 domains; ignoring "+str(len(seq_ws)-len(seqs_ws_4_domains))+" sequences that do not have 4 domains."))
    results = [global_thread_pool.apply_async(pfunc, args=( wc(s.split(' ')[2])+'TTTTT'+wc(s.split(' ')[1]),temp_in_C) ) for s in seqs_ws_4_domains]
    energies = [result.get() for result in results]    
    return energies
  

#############################################################################
# Run nupack's pfunc(), mfe(), binding(), and ViennaRNA's RNAduplex
#############################################################################

@lru_cache(maxsize=1000000)
def binding(seq1,seq2,temperature=53.0):
    """Computes the (partition function) free energy of association between two strands."""
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are 
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    # ddG_reaction == dG(products) - dG(reactants)
    return pfunc((seq1,seq2),temperature) - (pfunc(seq1,temperature) + pfunc(seq2,temperature))

@lru_cache(maxsize=1000000)
def duplex(seq,temp_in_C):
    """Computes the (partition function) free energy of a duplex relative to the single-stranded states."""
    return pfunc((seq, wc(seq)),temp_in_C) - pfunc((seq),temp_in_C) - pfunc((wc(seq)),temp_in_C)

def mfe_binding_async(seq1,seq2,temperature=53.0):
    """Computes the mfe of two strands."""
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are 
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    return mfe((seq1,seq2),temperature,unique_filename_param=seq1+'_'+seq2) 

def mfe_binding(seq1,seq2,temperature=53.0):
    """Computes the mfe of two strands."""
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are 
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    return mfe((seq1,seq2),temperature,unique_filename_param=seq1+'_'+seq2) 

def dGadjust(temperature,seqlen):
    R = 0.0019872041 # Boltzmann's constant in kcal/mol/K
    water_conc = 55.14 # molar concentration of water at 37 C; ignore temperature dependence, ~5%
    K = temperature + 273.15 # Kelvin
    adjust = R*K*math.log(water_conc) # converts from NUPACK mole fraction units to molar units, per association
    return adjust*(seqlen-1)

@lru_cache(maxsize=1000000)
def pfunc(seqtuple,temperature,adjust=True):
  """Calls NUPACK's pfunc on a complex consisting of the unique strands in
  seqtuple, returns dG.  temperature is in Celsius. More negative result (energy) is more favourable."""
  
  if type(seqtuple) is str:
    seqtuple = (seqtuple,)
  user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + ' '.join(map(str,list(range(1,len(seqtuple)+1))))
  p=sub.Popen(['pfunc','-T',str(temperature),'-multi','-material','dna'],
               stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,
               encoding='utf8')  # needed for python3 (as opposed to python2)             
  try:
    output = p.communicate(user_input)[0]
  except BaseException as error:
    p.kill()
    raise error
  
  lines = output.split('\n')
  
  if lines[-4] != "% Free energy (kcal/mol) and partition function:" :
    raise NameError('NUPACK output parsing problem')
  
  dG_str = lines[-3].strip()
  if dG_str.lower() == 'inf':
    # this occurs when two strands have MFE completely unpaired; should be 0 energy
    dG = 0.0
  else:
    if MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE: sign=1.0
    else: sign=-1.0 
    dG = sign*float(dG_str)
  
  if adjust:
    dG += dGadjust(temperature,len(seqtuple))    
  return dG


@lru_cache(maxsize=1000000)
def mfe(seqtuple,temperature,adjust=True,unique_filename_param=''):
  """Calls NUPACK's mfe on a complex consisting of the unique strands in
  seqtuple, returns dG.  temperature is in Celsius."""
      
  if type(seqtuple) is str:
    seqtuple = (seqtuple,)   
  file_data = str(len(seqtuple)) + '\n'
  for seq in seqtuple:
    file_data += seq + '\n'
  for i in range(len(seqtuple)):
    file_data += str(i+1) + ' '

  fname = '.'+unique_filename_param+'mfe_tmp_file'
  
  with open(fname+'.in','w') as f:
    f.write(file_data) 
    f.close()
  # print os.path.dirname(os.path.realpath(__file__))
  
  user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + ' '.join(map(str,list(range(1,len(seqtuple)+1))))
  p=sub.Popen(['mfe','-T',str(temperature),'-multi','-material','dna', fname],
               stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,
               encoding='utf8')  # added for python3             
  try:
    output = p.communicate(user_input)[0]
  except BaseException as error:
    p.kill()
    raise error
  
  with open(fname+'.mfe') as f:
    lines = f.readlines()
  
  if lines[1].strip() != "% Program: mfe" :
    raise NameError('NUPACK output parsing problem')
  
  if len(lines) == 12: # When two strands have MFE completely unpaired
    dG = 0.0
  else:
    if MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE: sign=1.0
    else: sign=-1.0 
    dG = sign*float(lines[14].strip())

  os.remove(fname+'.in') 
  os.remove(fname+'.mfe') 

  if adjust:
    dG += dGadjust(temperature,len(seqtuple))    
  return dG




def RNAduplex_multiple(seqpairs, temperature_in_C=53.0,NA_parameter_set=''): 
    """Calls RNAduplex on a list of pairs, specifically:
    [ (seq1, seq2), (seq2, seq3), (seq4, seq5), ... ]
    where seqi is a string over {A,C,T,G}. Temperature is in Celsius.
    Returns a list (in the same order as seqpairs) of negation of free energy 
    (so that more favourable means more positive)."""
    
    # NB: the string parameter_set needs to be exactly the intended filename; 
    # e.g. any extra whitespace characters causes RNAduplex to default to RNA parameter set without warning the user!

    if NA_parameter_set=='':
      NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 viennaRNA_PARAMETER_SET_DIRECTORY+'dna_mathews1999.par')    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
    else:
      NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 viennaRNA_PARAMETER_SET_DIRECTORY+NA_parameter_set)    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric

    # parameter_set = 'dna_mathews2004.par'  # Loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
    # parameter_set = 'dna_mathews1999.par'    # Gives better agreement with nupack than dna_mathews2004.par
    if not os.path.isfile(NA_parameter_set):
      raise ValueError('RNAduplex error: Error reading parameter file: ' + NA_parameter_set)
    
    # process the input into a string
    user_input = '\n'.join(seqpair[0]+'\n'+seqpair[1] for seqpair in seqpairs) + '@'

    p=sub.Popen(['RNAduplex','-P', NA_parameter_set, '-T', str(temperature_in_C), '--noGU'], 
                 stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,
                 encoding='utf8')  # added for python3
    
    try: output, stderr = p.communicate(user_input)
    except BaseException as error:
      p.kill()
      raise error    
    if stderr != '': # parsing error from RNAduplex 
      print(stderr) 
      if stderr.split('\n')[0] != 'WARNING: stacking enthalpies not symmetric':
        print('Stopping RNAduplex from loading default (RNA) parameter set')
        raise ValueError('RNAduplex error: Error reading parameter file ' + NA_parameter_set)  

    lines = output.split('\n')
    dG_list = []
    for line in lines[:-1]:
      if MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE: sign=1.0
      else: sign=-1.0 
      #dG_list.append(-float(line.split(':')[1].split('(')[1].split(')')[0]))
      dG_list.append(sign*float(line.split(':')[1].split('(')[1].split(')')[0]))
    return dG_list  # returns negated energies (i.e. more positive is more favourable)

def RNAduplex(seq1, seq2, temperature_in_C=53.0,NA_parameter_set=''): 
    """Calls RNAduplex on a pair of sequences. temperature is in Celsius.
    Returns negation of normal free energy so that result will be positive."""

    if NA_parameter_set=='':
      NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 viennaRNA_PARAMETER_SET_DIRECTORY+'dna_mathews1999.par')    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
    else:
      NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 viennaRNA_PARAMETER_SET_DIRECTORY+NA_parameter_set)    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
      
    user_input = seq1 +'\n' + seq2 + '@'
    if not os.path.isfile(NA_parameter_set):
        raise ValueError('RNAduplex error: Error reading parameter file ' + NA_parameter_set)
    
    # NB: the string parameter_set needs to be exactly the intended filename; 
    # any extra whitespace characters causes RNAduplex to default to RNA parameter set
    # which will raise an error
    #parameter_set = 'dna_mathews2004.par'
    #parameter_set = 'dna_mathews1999.par'
    
    p=sub.Popen(['RNAduplex','-P', NA_parameter_set, '-T', str(temperature_in_C), '--noGU'], 
                 stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE,
                 encoding='utf8')  # added for python3
    
    try: output, stderr = p.communicate(user_input)
    except BaseException as error:
        p.kill()
        raise error
    
    if stderr != '': # parsing error from RNAduplex 
        print(stderr) 
        if stderr.split('\n')[0] != 'WARNING: stacking enthalpies not symmetric':
            print('Stopping RNAduplex from loading default (RNA) parameter set')
            raise ValueError('RNAduplex error: Error reading parameter file ' + parameter_set )  
    dG_string = output.split(':')[1].split('(')[1].split(')')[0]
    if MORE_NEGAIVE_ENERGY_IS_MORE_FAVOURABLE: sign=1.0
    else: sign=-1.0 
    return sign*float(dG_string)

def hairpin(seq,temperature):
    """Computes the (partition function) free energy of single-strand secondary structure."""
    return pfunc((seq,),temperature)

def eval_colocated_end_pair(end1, end2, temperature):
    '''function used in the dsign code to compute binding energy for
    two co-located domains binding to each other.'''
    return max(hairpin(end1+end2, temperature), 
               min(hairpin(end1+'T'*4+end2, temperature), 
                   hairpin(end1+'A'*4+end2, temperature))) 
#    return max(hairpin(end1+end2, temperature), 
#               hairpin(end1+'T'*4+end2, temperature)) 


def idt_format_line(line):
  print(line)
  return len(line.split(","))==4 and set(line.split(",")[1]).issubset(set("ATCG atcg /iBiodT/"))

def seq_designer_format_line(line):
  return len(line.split("  "))==2 and set(line.split("  ")[1]).issubset(set("ATCG atcg /iBiodT/"))

def determine_file_type(f):
  lines = get_lines_from_file(f)   # removes "#"-comments
  consistent_with_idt_file = 1
  consistent_with_seq_designer_format_file = 1
  for line in lines:
    if not ( idt_format_line(line) ): consistent_with_idt_file=0
    if not ( seq_designer_format_line(line) ): consistent_with_seq_designer_format_file=0
  
  if consistent_with_seq_designer_format_file: 
    print("Input file is in format of sequence designer.")
    return "Sequence designer format"
  elif consistent_with_idt_file: 
    print("Input file is in IDT format.")
    return "IDT format"
  else: 
    print("Input file is in unknown format.")
    return "unknown format"

def sequencer_designer_file_format_to_IDT_file_format(input_filename):
  try:
    in_file = open(input_filename, "r")
  except IOError:
    print("file does not exist: " + input_filename); exit()
  print("processing file: "+input_filename)

  output_filename = os.path.splitext(input_filename)[0]+".idt"
  out_file = open(output_filename, "w")

  for line in in_file:
    l=line.strip()

    if l.startswith("#"):
      out_file.write(l+"\n") 

    if not l.startswith("#") and not l=="":
      name, sequence = l.split("  ")


      new_line = name+","+sequence
      if "/iBiodT/" in sequence: new_line+=",100nm,HPLC"
      else: new_line+=",25nm,STD"
      #print(new_line)
      out_file.write(new_line+"\n") 

  if "/iBiodT/" in ''.join(line for line in open(input_filename, "r")):
    print("Some strands have biotins (/iBiodT/), those were placed at 100nm scale in the output IDT-format file (all other strands at 25nm scale)")
  print("\nGenerated IDT-format file with filename: "+output_filename+"\n")
  return output_filename



def main(files,temp_in_C,domain_analysis=1,strand_analysis=1,tile_pair_check=1,
  lattice_binding_analysis=1,mfe_analysis=1,non_standard_tiles=0,pickle_data=0,num_bins=-1):
  # files is a list of paths to files with sequences in the appropriate format
  for f in files: 
    # Here, the following two lines are executed only for the purpose of syntax checking the input files.
    
    input_format = determine_file_type(f)
    if input_format == "Sequence designer format":
      #exit("TODO: Convert file format and write IDT file")
      f = sequencer_designer_file_format_to_IDT_file_format(f)
    elif input_format == "unknown format":
      exit("Error: Unknown file format")

    domains, domains_with_biotin = read_domains_from_idt_order(f,non_standard_tiles)
    names, seqs, seqs_bt, seq_ws, seq_ws_bt = read_strand_data_from_idt_order(f)
    
    if len(domains) != len(domains_with_biotin): 
      exit('  Parsing error: there are '+str(len(domains))+ ' domains when biotins are stripped out, and '+str(len(domains_with_biotin))+' when biotins not stripped out!')

    _analysis_directory = f +'__analysis/'+str(temp_in_C)+' C/'
    if not os.path.exists(_analysis_directory):
      os.makedirs(_analysis_directory)

    run_analysis(f,directory=_analysis_directory,temp_in_C=temp_in_C,domain_analysis=domain_analysis,
      strand_analysis=strand_analysis,tile_pair_check=tile_pair_check,
      lattice_binding_analysis=lattice_binding_analysis,run_mfe=mfe_analysis,
      non_standard_tiles=non_standard_tiles,pickle_data=pickle_data,num_bins=num_bins)


if __name__ == "__main__":
  args = parse_args()
  files = args.files
  domain_analysis = not(type(args.no_domains)==bool and args.no_domains)
  tile_pair_check = not(type(args.no_tile_pairs)==bool and args.no_tile_pairs) 
  strand_analysis = not(type(args.no_strand_analysis)==bool and args.no_strand_analysis)
  is_algorithmic_tile_set = not(type(args.no_lattice_binding)==bool and args.no_lattice_binding)
  lattice_binding_analysis = is_algorithmic_tile_set 
  mfe_analysis = not(type(args.no_mfe_analysis)==bool and args.no_mfe_analysis) 
  non_standard_tiles = type(args.non_standard_tiles)==bool and args.non_standard_tiles 
  if args.temperature: temp_in_C = args.temperature
  else: temp_in_C = _default_temp_in_C
  pickle_data = type(args.pickle)==bool and args.pickle
  if args.bins: num_bins = int(args.bins)
  else: num_bins=-1

  main(files=files,temp_in_C=temp_in_C,domain_analysis=domain_analysis,
      strand_analysis=strand_analysis,tile_pair_check=tile_pair_check,
      lattice_binding_analysis=lattice_binding_analysis,mfe_analysis=mfe_analysis,
      non_standard_tiles=non_standard_tiles,pickle_data=pickle_data,num_bins=num_bins)


