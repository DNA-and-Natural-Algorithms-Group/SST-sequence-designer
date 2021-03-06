'''
Shipped with DNA single-stranded tile (SST) sequence designer used in the following publication.
 "Diverse and robust molecular algorithms using reprogrammable DNA self-assembly"
 Woods*, Doty*, Myhrvold, Hui, Zhou, Yin, Winfree. (*Joint first co-authors)
 Nature, 2019

Generally sst_dsd processes Python 'ACTG' strings (as opposed to numpy arrays which are processed by dsd).
'''

from __future__ import print_function

import numpy as np
import math, string, random, sys, os
import subprocess as sub
from multiprocessing.pool import ThreadPool
import itertools

# from lru_cache import lru_cache # this was used when we had a hand-rolled lru_cache, but it's now in Python 3
try:
    from functools import lru_cache 
except ImportError:
    # in case this is python2
    from lru_cache import lru_cache

global_thread_pool = ThreadPool()

#viennaRNA_PARAMETER_SET_DIRECTORY = 'share/ViennaRNA/'
try:  
   os.environ["VIENNARNA_PARAMS_PATH"]
except KeyError: 
   print("Please set the environment path variable VIENNARNA_PARAMS_PATH to the location of the ViennaRNA parameters (dna_mathews1999.par and dna_mathews2004.par), ")
   print("For example by typing the following command (or putting it in your shell startup routine, e.g. for bash users in .bashrc):")
   print("export VIENNARNA_PARAMS_PATH=/usr/local/share/ViennaRNA/")
   sys.exit(1)

viennaRNA_PARAMETER_SET_DIRECTORY = os.environ.get('VIENNARNA_PARAMS_PATH')
DEFAULT_viennaRNA_PARAMETER_SET = viennaRNA_PARAMETER_SET_DIRECTORY+'dna_mathews1999.par'

# unix path must be able to find NUPACK, and NUPACKHOME must be set, as described in NUPACK installation instructions.

def dGadjust(temperature,seqlen):
    R = 0.0019872041 # Boltzmann's constant in kcal/mol/K
    water_conc = 55.14 # molar concentration of water at 37 C; ignore temperature dependence, ~5%
    K = temperature + 273.15 # Kelvin
    adjust = R*K*math.log(water_conc) # converts from NUPACK mole fraction units to molar units, per association
    return adjust*(seqlen-1)

@lru_cache(maxsize=1000000)
def pfunc(seqtuple, temperature, negate=False):
    """Calls NUPACK's pfunc on a complex consisting of the unique strands in
    seqtuple, returns dG.  temperature is in Celsius.

    If negate==True, returns negation of normal free energy so that result will be positive."""
    if type(seqtuple) is str:
        seqtuple = (seqtuple,)
    user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + ' '.join(map(str,list(range(1,len(seqtuple)+1))))

    p=sub.Popen(['pfunc','-T',str(temperature),'-multi','-material','dna'],
                 stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)             

    try:
        #output = p.communicate(user_input)[0]
        output = p.communicate(user_input.encode())[0]
        output = output.decode()
    except BaseException as error:
        p.kill()
        raise error
    
    #lines = output.split('\n')
    lines = output.split('\n')

    
    if lines[-4] != "% Free energy (kcal/mol) and partition function:" :
        raise NameError('NUPACK output parsing problem')
    
    dG_str = lines[-3].strip()
    if dG_str.lower() == 'inf':
        # this can occur when two strands have MFE completely unpaired; should be 0 energy
        dG = 0.0
    else:
        dG = float(dG_str)
    
    dG += dGadjust(temperature,len(seqtuple))
    if negate:
        dG = -dG
    return dG

def pfunc_multiple(seqtuples, temperature):
    """Calls NUPACK's pfunc on a several complexes consisting of the unique strands in
    seqtuple, returns dG.  temperature is in Celsius.

    Returns negations of normal free energy so that result will be positive."""
    if type(seqtuples[0]) is str:
        seqtuples = [(seqtuple,) for seqtuple in seqtuples]
        
    user_input = '\n'.join( str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + ' '.join(map(str,list(range(1,len(seqtuple)+1)))) for seqtuple in seqtuples ) + '\n-1\n'
    
    p=sub.Popen(['pfunc_multi','-T',str(temperature),'-multi','-material','dna'],
                 stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)   
              
    try:
        output = p.communicate(user_input.encode())[0]
        output = output.decode()
    except BaseException as error:
        p.kill()
        raise error
        
    lines = output.split('\n')[6:-1:2]
        
    dGa = dGadjust(temperature, len(seqtuples[0]))
    pfunc_energies = [-(float(line)+dGa) for line in lines]

    return pfunc_energies


@lru_cache(maxsize=1000000)
def mfe(seqtuple, temperature, adjust=True, unique_filename_param='', negate=False):
    """Calls NUPACK's mfe on a complex consisting of the unique strands in
    seqtuple, returns dG.  temperature is in Celsius."""

    if type(seqtuple) is str:
        seqtuple = (seqtuple,)
    file_data = str(len(seqtuple)) + '\n'
    for seq in seqtuple:
        file_data += seq + '\n'
    for i in range(len(seqtuple)):
        file_data += str(i + 1) + ' '

    fname = '.' + unique_filename_param + 'mfe_tmp_file'

    with open(fname + '.in', 'w') as f:
        f.write(file_data)
        f.close()
    # print os.path.dirname(os.path.realpath(__file__))

    user_input = str(len(seqtuple)) + '\n' + '\n'.join(seqtuple) + '\n' + ' '.join(
        map(str, list(range(1, len(seqtuple) + 1))))
    p = sub.Popen(['mfe', '-T', str(temperature), '-multi', '-material', 'dna', fname],
                  stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE,
                  universal_newlines=True)

    try:
        output = p.communicate(user_input)[0]
    except BaseException as error:
        p.kill()
        raise error

    with open(fname + '.mfe') as f:
        lines = f.readlines()

    if lines[1].strip() != "% Program: mfe":
        raise NameError('NUPACK output parsing problem')

    if len(lines) == 12:  # When two strands have MFE completely unpaired
        dG = 0.0
    else:
        if not negate:
            sign = 1.0
        else:
            sign = -1.0
        dG = sign * float(lines[14].strip())

    os.remove(fname + '.in')
    os.remove(fname + '.mfe')

    if adjust:
        dG += dGadjust(temperature, len(seqtuple))
    return dG


def mfe_binding(seq1,seq2,temperature=53.0, negate=False):
    """Computes the mfe of two strands."""
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    return mfe((seq1,seq2),temperature,unique_filename_param=seq1+'_'+seq2, negate=negate)


def RNAduplex_multiple(seqpairs, temperature_in_C, NA_parameter_set, negate=False):
    """Calls RNAduplex on a list of pairs, specifically:
    [ (seq1, seq2), (seq2, seq3), (seq4, seq5), ... ]
    where seqi is a string over {A,C,T,G}. Temperature is in Celsius.
    Returns a list (in the same order as seqpairs) of negation of free energy 
    (so that more favourable means more positive)."""
    
    # NB: the string NA_parameter_set needs to be exactly the intended filename; 
    # e.g. any extra whitespace characters causes RNAduplex to default to RNA parameter set without warning the user!
    
    if NA_parameter_set=='':
        NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 DEFAULT_viennaRNA_PARAMETER_SET)    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
    else:
        NA_parameter_set = os.path.join(os.path.dirname(__file__),
                                 viennaRNA_PARAMETER_SET_DIRECTORY+NA_parameter_set)    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric

    if not os.path.isfile(NA_parameter_set):
        raise ValueError('RNAduplex error: Error reading parameter file: ' + NA_parameter_set)

    # process the input into a string
    user_input = '\n'.join(seqpair[0]+'\n'+seqpair[1] for seqpair in seqpairs) + '\n@\n'

    got_results = False
    dG_list = []
    while not got_results:
        # When porting the code from python2 to python3 we found an issue with sub.Popen(). 
        # Passing either of the keyword arguments universal_newlines=True or encoding='utf8' 
        # solves the problem for python3.6. For python3.7 (but not 3.6) one can use text=True
        p=sub.Popen(['RNAduplex','-P', NA_parameter_set, '-T', str(temperature_in_C), '--noGU'], # , '--noconv' make sense to use this, but it's untested by us 
                     stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)  

        try: 
          output, stderr = p.communicate(user_input.encode())
          output = output.decode()
          stderr = stderr.decode()
        except BaseException as error:
            p.kill()
            raise error
        dG_list = []
        if stderr.strip() != '': # parsing error from RNAduplex
            print('Warning or error from RNAduplex: ', stderr)
            if stderr.split('\n')[0] != 'WARNING: stacking enthalpies not symmetric':
                print('Stopping RNAduplex from loading default (RNA) parameter set')
                # print('Re-calling RNAduplex due to an error (using DNA parameter file ' + NA_parameter_set +')')
                print('RNAduplex says:' + str(stderr.split('\n')[0]))
                raise ValueError('RNAduplex error: Error reading parameter file: ' + NA_parameter_set)

        lines = output.split('\n')
        for line in lines[:-1]:
            energy = float( line.split(':')[1].split('(')[1].split(')')[0] )
            if negate:
                energy = -energy
            dG_list.append(energy)
        if len(lines) - 1 != len(seqpairs):
            raise ValueError('lengths do not match: #lines:{} #seqpairs:{}'.format(len(lines)-1, len(seqpairs)))
        got_results = True
    
    return dG_list  # returns negated energies (i.e. more positive is more favourable)

# def RNAcofold_multiple(seqpairs, temperature_in_C):
#     """Calls RNAduplex on a list of pairs, specifically:
#     [ (seq1, seq2), (seq2, seq3), (seq4, seq5), ... ]
#     where seqi is a string over {A,C,T,G}. Temperature is in Celsius.
#     Returns a list (in the same order as seqpairs) of negation of free energy
#     (so that more favourable means more positive)."""
#
#     # NB: the string parameter_set needs to be exactly the intended filename;
#     # e.g. any extra whitespace characters causes RNAduplex to default to RNA parameter set without warning the user!
#
#     # parameter_set = 'dna_mathews2004.par'  # Loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
#     #parameter_set = 'dna_mathews1999.par'    # Gives better agreement with nupack than dna_mathews2004.par
#     parameter_set = os.path.join(os.path.dirname(__file__),
#                                  'nupack_viennaRNA/dna_mathews1999.par')    # Gives better agreement with nupack than dna_mathews2004.par. Note that loading parameter set dna_mathews2004.par throws a warning encoded in that parameter set:  WARNING: stacking enthalpies not symmetric
#
#
#     # process the input into a string
#     user_input = '\n'.join(seqpair[0]+'&'+seqpair[1] for seqpair in seqpairs) + '\n@\n'
#
#     p=sub.Popen(['RNAcofold','-P', parameter_set, '-T', str(temperature_in_C), '--noGU', '--noconv', '-p'],
#                  stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
#
#     try:
#         output, stderr = p.communicate(user_input.encode())
#         output = output.decode()
#         stderr = stderr.decode()
#     except BaseException as error:
#         p.kill()
#         raise error
#     print('output="', output, '"')
#     print('stderr="', stderr, '"')
#     if stderr != '': # parsing error from RNAduplex
#         print(stderr)
#         if stderr.split('\n')[0] != 'WARNING: stacking enthalpies not symmetric':
#             print('Stopping RNAduplex from loading default (RNA) parameter set')
#             raise ValueError('RNAduplex error: Error reading parameter file ' + parameter_set)
#     return
#
#     #DD: the code below is unreachable; delete it (or delete this whole function, since it doesn't appear to be called in any file)
#     lines = output.split('\n')
#     dG_list = []
#     for line in lines[:-1]:
#         dG_list.append(-float( line.split(':')[1].split('(')[1].split(')')[0] ))
#     if len(lines) - 1 != len(seqpairs):
#         raise ValueError('lengths do not match: #lines:{} #seqpairs:{}'.format(len(lines)-1, len(seqpairs)))
#     return dG_list  # returns negated energies (i.e. more positive is more favourable)


# maketrans is in string if python2 but in str if python3
try:
    _wctable = str.maketrans('ACGTacgt','TGCAtgca')
except AttributeError:
    _wctable = string.maketrans('ACGTacgt','TGCAtgca')

def wc(seq):
    """Return reverse Watson-Crick complement of seq"""
    return seq.translate(_wctable)[::-1]

def duplex(seq, temperature, subtract_indv=True, negate=False):
    """Computes the (partition function) free energy of a duplex."""
    seq1 = seq
    seq2 = wc(seq)
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are 
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    association_energy = pfunc((seq1,seq2), temperature, negate)
    if subtract_indv:
        # ddG_reaction == dG(products) - dG(reactants)
        association_energy -= (pfunc(seq1, temperature, negate) + pfunc(seq2, temperature, negate))
    return association_energy

def hairpin(seq,temperature,negate=False):
    """Computes the (partition function) free energy of single-strand secondary structure."""
    return pfunc((seq,), temperature, negate)

def binding(seq1, seq2, temperature, negate=False):
    """Computes the (partition function) free energy of association between two strands."""
    # this is a hack to save time since (seq1,seq2) and (seq2,seq1) are 
    #   considered different tuples hence are cached differently by lrucache;
    #   but pfunc is a symmetric function so it's safe to swap the order
    if seq1 > seq2: seq1,seq2 = seq2,seq1
    # ddG_reaction == dG(products) - dG(reactants)
    return pfunc((seq1,seq2), temperature, negate) - (pfunc(seq1, temperature, negate) + pfunc(seq2, temperature, negate))

def randomseq(length,bases='ACTG'):
    """Chooses a random DNA sequence."""
    return ''.join(random.choice(bases) for _ in range(length))

def domain_equal_strength(seq, temperature, lowPF, highPF):
    """test roughly equal strength of domains (according to partition function)"""
    #dG = duplex(seq,temperature)
    negate = True
    dG = binding(seq, wc(seq), temperature, negate)
    return lowPF <= dG <= highPF

def domain_no_sec_struct(seq,temperature,individual,threaded):
    """test lack of secondary structure in domains"""
    negate = True
    if threaded:
        results = [global_thread_pool.apply_async(hairpin, args=(s,temperature)) for s in (seq, wc(seq))]
        e_seq,e_wcseq = [result.get() for result in results]
        return e_seq <= individual and e_wcseq <= individual
    else:
        return hairpin(seq, temperature, negate) <= individual and hairpin(wc(seq), temperature, negate) <= individual


LOG_ENERGY = False
def log_energy(energy):
    if LOG_ENERGY:
        print('%.1f' % energy)

def domain_orthogonal(seq,seqs,temperature,orthogonality,orthogonality_ave=-1,threaded=True):
    """test orthogonality of domain with all others and their wc complements"""
    negate = True
    if threaded:
        results = [global_thread_pool.apply_async(binding, args=(s,s,temperature)) for s in (seq, wc(seq))]
        energies = [result.get() for result in results]
        if max(energies) > orthogonality: 
            return False
    else:
        ss = binding(seq, seq, temperature, negate)
        log_energy(ss)
        if ss > orthogonality: 
            return False
        wsws = binding(wc(seq), wc(seq), temperature, negate)
        log_energy(wsws)
        if wsws > orthogonality: 
            return False
    energy_sum = 0.0
    for altseq in seqs:
        if threaded:
            results = [global_thread_pool.apply_async(binding, args=(seq1,seq2,temperature)) 
                       for seq1,seq2 in itertools.product((seq, wc(seq)), (altseq, wc(altseq)))]
            energies = [result.get() for result in results]
            if max(energies) > orthogonality: 
                return False
            energy_sum += sum(energies)
        else:
            sa = binding(seq, altseq, temperature, negate)
            log_energy(sa)
            if sa > orthogonality: 
                return False
            sw = binding(seq, wc(altseq), temperature, negate)
            log_energy(sw)
            if sw > orthogonality: 
                return False
            wa = binding(wc(seq), altseq, temperature, negate)
            log_energy(wa)
            if wa > orthogonality: 
                return False
            ww = binding(wc(seq), wc(altseq), temperature, negate)
            log_energy(ww)
            if ww > orthogonality: 
                return False
            energy_sum += sa+sw+wa+ww
    if orthogonality_ave > 0:
        energy_ave = energy_sum / (4*len(seqs)) if len(seqs) > 0 else 0.0
        return energy_ave <= orthogonality_ave
    else:
        return True

def domain_pairwise_concatenated_no_sec_struct(seq,seqs,temperature,concat,concat_ave=-1,threaded=True):
    """test lack of secondary structure in concatenated domains"""
    energy_sum = 0.0
    negate = True
    for altseq in seqs:
        wc_seq = wc(seq)
        wc_altseq = wc(altseq)
        if threaded:
            results = [global_thread_pool.apply_async(hairpin, args=(seq1+seq2, temperature)) for (seq1,seq2) in
                       [(seq,altseq), 
                        (seq,wc_altseq), 
                        (wc_seq,altseq), 
                        (wc_seq,wc_altseq), 
                        (altseq,seq), 
                        (wc_altseq,seq), 
                        (altseq,wc_seq), 
                        (wc_altseq,wc_seq)]]
            energies = [result.get() for result in results]
#             print len(results)
#             print 'pair: %s' % [round(e,1) for e in energies]
            if max(energies) > concat: return False
            energy_sum += sum(energies)
        else:
            seq_alt = hairpin(seq + altseq, temperature, negate)
            if seq_alt > concat: return False
            seq_wcalt = hairpin(seq + wc_altseq, temperature, negate)
            if seq_wcalt > concat: return False
            wcseq_alt = hairpin(wc_seq + altseq, temperature, negate)
            if wcseq_alt > concat: return False
            wcseq_wcalt = hairpin(wc_seq + wc_altseq, temperature, negate)
            if wcseq_wcalt > concat: return False
            alt_seq = hairpin(altseq + seq, temperature, negate)
            if alt_seq > concat: return False
            alt_wcseq = hairpin(altseq + wc_seq, temperature, negate)
            if alt_wcseq > concat: return False
            wcalt_seq = hairpin(wc_altseq + seq, temperature, negate)
            if wcalt_seq > concat: return False
            wcalt_wcseq = hairpin(wc_altseq + wc_seq, temperature, negate)
            if wcalt_wcseq > concat: return False
            energy_sum += (seq_alt + seq_wcalt + wcseq_alt + wcseq_wcalt + 
                           alt_seq + alt_wcseq + wcalt_seq + wcalt_wcseq)
    if concat_ave > 0:
        energy_ave = energy_sum / (8*len(seqs)) if len(seqs) > 0 else 0.0
        return energy_ave <= concat_ave
    else:
        return True

def check(s1, s2, T, orthogonality, pass_tests, idx):
    negate = True
    pass_tests[idx] = (binding(s1, s2, T, negate) <= orthogonality)

# maketrans is in string if python2 but in str if python3
try:
    _binaryGCTable = str.maketrans('ACTG','0101')
except AttributeError:
    _binaryGCTable = string.maketrans('ACTG','0101')

def domain_concatenated_no4GC(seq,seqs):
    """prevent {G,C}^4 under concatenation"""
    for altseq in seqs:
        catseq = altseq+seq+altseq
        strength = catseq.translate(_binaryGCTable)
        if '1111' in strength: return False
    return True

def domain_no4GC(seq):
    """prevent {G,C}^4"""
    return '1111' not in seq.translate(_binaryGCTable)

def domain_concatenated_no4Gor4C(seq,seqs):
    """prevent G^4 and C^4 under concatenation"""
    for altseq in seqs:
        catseq = altseq+seq+altseq
        if 'GGGG' in catseq: 
#             print '|GGGG# seq: %s altseq: %s|' % (seq,altseq)
            return False
        if 'CCCC' in catseq: 
#             print '|CCCC# seq: %s altseq: %s|' % (seq,altseq)
            return False
    return True

def all_cats(seqarr,seqsarr):
    """Return all sequences obtained by concatenating seqarr to either end of
    a sequence in seqsarr.

    For example, all_cats([0,1,2,3], [[3,3,3],
                                      [0,0,0]]) returns
    [[0,1,2,3,3,3,3],
    [3,3,3,0,1,2,3],
    [0,1,2,3,0,0,0],
    [0,0,0,0,1,2,3]]"""
    seqarr=np.asarray([seqarr])
    seqsarr=np.asarray(seqsarr)
    ar=seqarr.repeat(seqsarr.shape[0],axis=0)
    ret=np.concatenate((seqsarr,ar),axis=1)
    ret2 = np.concatenate((ar,seqsarr),axis=1)
    ret = np.concatenate((ret,ret2))
    return ret

def log_bad_end(reason, log):
    if log:
        sys.stdout.write(reason)
        sys.stdout.flush()

# def has_hairpin(seq,hairpin):
#     for i in range(len(seq) - 2*hairpin - 3):
#         sub = seq[i:i+hairpin]
#         subWC = wc(sub)
#         if subWC in seq[i+hairpin+3]:
#             return False
#
# def domain_concatenated_no_hairpin(seq,seqs,hairpin=5):
#     """prevent hairpins of stem length 5 or more"""
#     for altseq in seqs:
#         catseq = altseq+seq
#         if has_hairpin(catseq, hairpin): return False
#         catseq = seq+altseq
#         if has_hairpin(catseq, hairpin): return False
#     return True
#
# def nextseq(init_seqs,new_seqs,iterator,temperature,lowPF,highPF,individual,
#             orthogonality,concat,orthogonality_ave,concat_ave,
#             prevent4GC,prevent4G4C,hairpin,threaded=True):
#     """Return next sequence from iterator that "gets along" with the sequences
#     already in existing_seqs according to parameters."""
#     all_seqs = init_seqs + new_seqs
#     log = True
#     num_searched = 0
#     sys.stdout.write('.')
#     sys.stdout.flush()
#     #seqsarr = dsd.seqs2arr(existing_seqs)
#     for seq in iterator:
#         num_searched += 1
# #         sys.stdout.write('.')
# #         sys.stdout.flush()
#         if wc(seq) in all_seqs: continue
#         if prevent4GC and not domain_no4GC(seq): # domain_concatenated_no4GC(seq,new_seqs):
#             log_bad_end('gc4_',log)
#             continue
#         if not prevent4GC and (prevent4G4C and ('GGGG' in seq or 'CCCC' in seq)): # domain_concatenated_no4Gor4C(seq,new_seqs)):
#             log_bad_end('g4c4_',log)
#             continue
#         #if not domain_concatenated_no_hairpin(seq,existing_seqs,hairpin): continue
#         #if hairpin and not domain_concatenated_no_hairpin_arr(seq,seqsarr,hairpin): continue
#         if not domain_equal_strength(seq,temperature,lowPF,highPF):
#             log_bad_end('eq_',log)
#             continue
#         if not domain_no_sec_struct(seq,temperature,individual,threaded):
#             log_bad_end('idv_', log)
#             continue
#         if not domain_pairwise_concatenated_no_sec_struct(seq,new_seqs,temperature,concat,concat_ave,threaded):
#             log_bad_end('cat_',log)
#             continue
#         if not domain_orthogonal(seq,all_seqs,temperature,orthogonality,orthogonality_ave,threaded):
#             log_bad_end('orth_',log)
#             continue
#         sys.stdout.write('.\n')
#         sys.stdout.flush()
#         return seq,num_searched
#     raise ValueError('no more sequences to search')
#
# def learnSL(lengths,lowPF,highPF,temperature,num_samples = 100):
#     """Learn appropriate upper and lower bounds for SantaLucia energy (as
#     calculated by DNASeqList.wcenergies) that preserve "many" sequences whose
#     binding energies according to binding function remain.
#
#     Current algorithm gets sample min and max and chooses lowSL and highSL endpoints
#     to contain the middle two quartiles.
#
#     Originally used minimum variance unbiased estimators for min and max,
#     example 2 here:
#     http://www.d.umn.edu/math/Technical%20Reports/Technical%20Reports%202007-/TR%202010/TR_2010_8.pdf
#
#     but that gave too large a range."""
#     inrange = 0
#     negate = True
#     print('Searching for optimal SantaLucia energy range within binding energy ' \
#           + 'lowSL %.2f and highSL %.2f\n***********************' % (lowPF,highPF))
#     energies = []
#     while inrange < num_samples:
#         sys.stdout.write('.')
#         sys.stdout.flush()
#         seq = randomseq(random.choice(lengths))
#         energyBinding = binding(seq, wc(seq), temperature, negate)
#         if lowPF <= energyBinding <= highPF:
#             inrange += 1
#             energySL = dsd.wcenergy(seq, temperature)
#             energies.append(energySL)
#             #print energySL
#         sys.stdout.write('.')
#         sys.stdout.flush()
#     print()
#     energies.sort()
#     #print [round(e,2) for e in energies]
#     lowerPos = num_samples//4
#     upperPos = num_samples - lowerPos
#     assert lowerPos < upperPos
#     lowSL = energies[lowerPos]
#     highSL = energies[upperPos]
#     return (lowSL,highSL)
#
# def learnPF(seqlists,temperature,num_samples = 100):
#     """Learn appropriate upper and lower bounds for partition function energy (as
#     calculated by sst_dsd.binding(s,wc(s))) that preserve "many" sequences whose
#     SantaLucia binding energies according to dsd.wcenergy function remain.
#
#     Current algorithm gets sample min and max and chooses low and high endpoints
#     to contain the middle two quartiles.
#
#     Originally used minimum variance unbiased estimators for min and max,
#     example 2 here:
#     http://www.d.umn.edu/math/Technical%20Reports/Technical%20Reports%202007-/TR%202010/TR_2010_8.pdf
#
#     but that gave too large a range."""
#     negate = True
#     energies = []
#     for i in range(num_samples):
#         seqlist = random.choice(seqlists)
#         if not (isinstance(seqlist,list) or isinstance(seqlist,dsd.DNASeqList)):
#             raise TypeError('seqlist must be DNASeqList or list, not %s' % seqlist.__class__)
#
#         # elements of seqlists could be either DNASeq objects or lists of strings
#         idx = random.randint(0,len(seqlist)-1)
#         seq = seqlist[idx]
#         energyPF = duplex(seq, temperature, negate)
#         energies.append(energyPF)
#         #print 'energyPF:%.2f energySL:%.2f' % (energyPF,seqlist.wcenergy(idx))
#         #sys.stdout.write('.')
#         #sys.stdout.flush()
#     energies.sort()
#     #print [round(e,2) for e in energies]
#     lowerPos = num_samples//6
#     upperPos = num_samples - lowerPos
#     assert lowerPos < upperPos
#     highPF = energies[upperPos]
#     lowPF = energies[lowerPos]
#     return (lowPF,highPF)
#
# def prefilter_length_10_11_autofilter(temperature):
#     raise NotImplementedError()
#     seqs10 = dsd.DNASeqList(length=10)
#     seqs11 = dsd.DNASeqList(length=11)
#     seqs10.sortSeqsByWCEnergy()
#     seqs11.sortSeqsByWCEnergy()
#     seqs10,seqs11 = dsd.filterByCommonDensityVectorized(seqs10, seqs11, delta=0.01)
#     seqs10 = seqs10.filterSeqsByGQuadCQuad()
#     seqs11 = seqs11.filterSeqsByGQuadCQuad()
#     return (seqs10.toList(), seqs11.toList())
#
# def prefilter_length_10_11(lowDG, highDG, temperature, endGC, convert_to_list=True):
#     """Return sequences of length 10 and 11 with wc energies between given values."""
#     s10 = dsd.DNASeqList(length=10)
#     s11 = dsd.DNASeqList(length=11)
#     s10 = s10.filter_energy(low=lowDG, high=highDG, temperature=temperature)
#     s11 = s11.filter_energy(low=lowDG, high=highDG, temperature=temperature)
#     #s10 = s10.filterSeqsByGQuadCQuad()
#     #s11 = s11.filterSeqsByGQuadCQuad()
#     forbidden_subs = ['%s%s%s%s' % (a,b,c,d) for a in ['G','C']
#                                              for b in ['G','C']
#                                              for c in ['G','C']
#                                              for d in ['G','C']]
#     s10 = s10.filter_substring(forbidden_subs)
#     s11 = s11.filter_substring(forbidden_subs)
#     if endGC:
#         print('Removing any domains that end in either A or T; also ensuring every domain has an A or T within 2 indexes of the end')
#         s10 = s10.filter_endGC()
#         s11 = s11.filter_endGC()
#     #s10.sortSeqsByWCEnergy()
#     #s11.sortSeqsByWCEnergy()
#     if convert_to_list:
#         s10 = s10.toList()
#         s11 = s11.toList()
#     for seqs in (s10,s11):
#         if len(seqs) == 0:
#             raise ValueError('lowDG %.2f and highDG %.2f too strict! no sequences of length %d found' % (lowDG,highDG,seqs.seqlen))
#     return (s10, s11)
#
# def design_domains_10_11(howmany=1, temperature=53.0, lowSL=None, highSL=None,
#                          lowPF=None, highPF=None, individual=1.0,
#                          orthogonality=4.5, concat=3.3,
#                          orthogonality_ave=4.5, concat_ave=3.3,
#                          prevent4GC=False, prevent4G4C=True,
#                          hairpin=0, pr=None, init_seqs=[], endGC=False):
#     """Like design domains but specialized to length 10 and 11 domains.
#     Also iterator uses custom code to start with a small(ish) set of sequences with
#     similar binding energies.
#
#     lowSL and highSL are lower and upper limits on the energy as reported by
#     the SantaLucia nearest neighbor energy model as computed in
#     DNASeqList.wcenergy(idx) (i.e., the energy of the sequence bound to its
#     complement). It should be related to target and spread,
#     which are energies related to the partition function of the ordered complex
#     consisting of the sequence and its complement, but frankly I'm not sure
#     what the relationship is in general. The SL energy is MUCH faster to
#     calculate on all sequences at once, so serves as a fast preliminary filter.
#     """
#     if (lowSL == None or highSL == None) and (lowPF == None or highPF == None):
#         raise ValueError('At least one of the pairs (lowSL,highSL) or (lowPF,highPF) must be specified.')
#     if lowSL:
#         if not highSL: raise ValueError('lowSL specified but not highSL')
#         print('Using user-specified SantaLucia energy range [%.2f,%.2f]' % (lowSL,highSL))
#     if lowPF:
#         if not highPF: raise ValueError('lowPF specified but not highPF')
#         print('Using user-specified partition energy range [%.2f,%.2f]' % (lowPF,highPF))
#
#     if lowSL == None or highSL == None:
#         print('learning SantaLucia energy')
#         lowSL,spreadSL_ret = learnSL((10,11),lowPF,highPF,temperature)
#         if not highSL: highSL = spreadSL_ret
#         print('Using learned SantaLucia energy range [%.2f,%.2f]' % (lowSL,highSL))
#         seqs10,seqs11 = prefilter_length_10_11(lowSL,highSL,temperature,endGC)
#     elif lowPF == None or highPF == None:
#         s10,s11 = prefilter_length_10_11(lowSL,highSL,temperature,endGC,convert_to_list=False)
#         print('learning partition energy')
#         lowPF,highPF = learnPF((s10,s11),temperature,num_samples=100)
#         print('Using learned partition energy range [%.2f,%.2f]' % (lowPF,highPF))
#         seqs10 = s10.toList()
#         seqs11 = s11.toList()
#     elif lowPF and lowSL and highPF and highSL:
#         seqs10,seqs11 = prefilter_length_10_11(lowSL,highSL,temperature,endGC)
#
#
#     print('num length-10 seqs found:%d' % len(seqs10))
#     print('num length-11 seqs found:%d' % len(seqs11))
#     random.shuffle(seqs10)
#     random.shuffle(seqs11)
#     new_seqs = []
#     on10 = True
#     iter10 = iter(seqs10)
#     iter11 = iter(seqs11)
#
#     num_total10 = len(seqs10)
#     num_total11 = len(seqs11)
#     num_searched10 = 0
#     num_searched11 = 0
#
#     if pr: pr.enable()
#     raise NotImplementedError('put description of constraint parameters in file name')
#     filename = 'seqs/sequences_%s.txt' % str(datetime.datetime.now()).replace(' ','_').replace(':','-')
#     with open(filename,'w') as f:
#         while len(new_seqs) < howmany:
#             iterator = iter10 if on10 else iter11
#             sys.stdout.write(str(len(new_seqs)))
#             start_time = time.time()
#             seq,num_searched = nextseq(init_seqs=init_seqs,new_seqs=new_seqs,iterator=iterator,
#                                temperature=temperature,lowPF=lowPF,highPF=highPF,
#                                individual=individual,
#                                orthogonality=orthogonality,concat=concat,
#                                orthogonality_ave=orthogonality_ave,concat_ave=concat_ave,
#                                prevent4GC=prevent4GC,prevent4G4C=prevent4G4C,hairpin=hairpin)
#             tot_time = time.time() - start_time
#             if on10: num_searched10 += num_searched
#             else: num_searched11 += num_searched
#             if seq:
#                 new_seqs.append(seq)
#                 print(seq)
#                 print((' time: %5.1f secs' % tot_time))
#                 print((' length 10 searched:  %6d' % num_searched10), end=' ')
#                 print((' length 10 remaining: %6d' % (num_total10-num_searched10)), end=' ')
#                 print((' length 11 searched:  %6d' % num_searched11), end=' ')
#                 print((' length 11 remaining: %6d' % (num_total11-num_searched11)))
#                 f.write(seq+'\n')
#                 f.flush()
#             else:
#                 print('Could not find %d sequences matching your criteria' % howmany)
#                 break
#             on10 = not on10
#     if pr: pr.disable()
#     return ([s for s in new_seqs if len(s) == 10],
#             [s for s in new_seqs if len(s) == 11])
#
# def main():
# #     init_seqs = ['TTGAGGAGAG',
# #                  'TGTAGTAGGC',
# #                  'ATGTTTTGGG',
# #                  'TTGGTGATTC',
# #                  'AGTTTGTTGC',
# #                  'ATAGTGGGAG',
# #                  'AAGGATGGAC',
# #                  'TGTAATTGGC',
# #                  'ATAGGGATGC',
# #                  'TGAGGGTTAG',
# #                  'TGAAATGGTC',
# #                  'TAAGTGGTGG',
# #                  'TGATGAGGTG',
# #                  'TGTGTGAGAC',
# #                  'AAGAGAGGAC',
# #                  'AGGATTGGAG']
#     init_seqs=['GCCTACTACA', 'CCATCAAACCA', 'GTCCTACACTT', 'GTCCTCTCTT', 'AAGTGTAGGAC',
#                'AAGAGAGGAC', 'AGGATTGGAG', 'AGAGATTGTTC', 'CTCCAATCCT', 'GAACAATCTCT',
#                'GCTACACAATT', 'CACCTCATCA', 'AATTGTGTAGC', 'TGATGAGGTG', 'TGTGTGAGAC',
#                'TTGAAGAAGAC', 'GTCTCACACA', 'GTCTTCTTCAA', 'CCAACCTATTT', 'GACCATTTCA',
#                'AAATAGGTTGG', 'TGAAATGGTC', 'TAAGTGGTGG', 'AGTAAGAAGGC', 'CCACCACTTA',
#                'GCCTTCTTACT', 'CTCACTACATT', 'GCATCCCTAT', 'AATGTAGTGAG', 'ATAGGGATGC',
#                'TGAGGGTTAG', 'TGGTAAGGAAC', 'CTAACCCTCA', 'GTTCCTTACCA', 'GCTCTTCACAA',
#                'GTCCATCCTT', 'TTGTGAAGAGC', 'AAGGATGGAC', 'TGTAATTGGC', 'TGGGATAGTAG',
#                'GCCAATTACA', 'CTACTATCCCA', 'GACTTATCCAA', 'GCAACAAACT', 'TTGGATAAGTC',
#                'AGTTTGTTGC', 'ATAGTGGGAG', 'AATTAGGTAGC', 'CTCCCACTAT', 'GCTACCTAATT',
#                'GCAATATCACA', 'CCCAAAACAT', 'TGTGATATTGC', 'ATGTTTTGGG', 'TTGGTGATTC',
#                'TATTGTTAGGC', 'GAATCACCAA', 'GCCTAACAATA', 'CCCTACAACAA', 'CTCTCCTCAA',
#                'TTGTTGTAGGG', 'TTGAGGAGAG', 'TGTAGTAGGC', 'TGGTTTGATGG', 'CTCCAATCCT',
#                'GAACAATCTCT', 'CCTCAAATACA', 'CCATCATCAA', 'TGTATTTGAGG', 'TTGATGATGG',
#                'TGTGTGAGAC', 'TTGAAGAAGAC', 'CCACCACTTA', 'GCCTTCTTACT', 'GACCTACCATA',
#                'CCTCAACTCA', 'TATGGTAGGTC', 'TGAGTTGAGG', 'TGAGGGTTAG', 'TGGTAAGGAAC',
#                'GCCAATTACA', 'CTACTATCCCA', 'GACAACTACCT', 'GCTCAATACA', 'AGGTAGTTGTC',
#                'TGTATTGAGC', 'ATAGTGGGAG', 'AATTAGGTAGC', 'GAATCACCAA', 'GCCTAACAATA',
#                'CACTAATCACA', 'CTCTCTACCA', 'TGTGATTAGTG', 'TGGTAGAGAG', 'TGTAGTAGGC',
#                'TGGTTTGATGG']
# #     init_seqs=[]
#     s10,s11 = design_domains_10_11(howmany=40,temperature=53.0,
#                                    lowSL=10.0,highSL=10.5,
#                                    #targetPF=9.0,spreadPF=0.5,
#                                    prevent4GC=False, prevent4G4C=True,
#                                    individual=1.0,
#                                    orthogonality=4.0,concat=2.5,
#                                    orthogonality_ave=2.2,concat_ave=1.2,
#                                    init_seqs=init_seqs, endGC=True)
#
#     delim = '*'*79
#     print(delim)
#     print('Python representation:')
#     print(delim)
#     print('s10 = %s' % s10)
#     print('s11 = %s' % s11)
#     print(delim)
#     print('Sequences delimited by newlines:')
#     print(delim)
#     for s in s10: print(s)
#     for s in s11: print(s)
#
# if __name__ == "__main__":
#     main()


