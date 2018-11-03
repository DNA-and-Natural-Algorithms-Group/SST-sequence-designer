'''
Library for doing sequence design that can be expressed as linear algebra
operations for rapid processing by numpy (e.g., generating all DNA sequences
of a certain length and calculating all their full duplex binding energies
in the nearest neighbor model and filtering those outside a given range).
'''

from __future__ import division
import numpy as np
import sst_dsd
from lru_cache import lru_cache
import itertools as it

PYPY = True

_SCIPY_AVAIL = not PYPY

_bits2base = ['A','C','G','T']
_base2bits = {'A':0b00, 'C':0b01, 'G':0b10, 'T':0b11,
              'a':0b00, 'c':0b01, 'g':0b10, 't':0b11}

def idx2seq(idx,length):
    '''Return the lexicographic idx'th DNA sequence of given length.'''
    seq = ['x']*length
    for i in range(length-1,-1,-1):
        seq[i] = _bits2base[idx & 0b11]
        idx >>= 2
    return ''.join(seq)

def seq2arr(seq):
    '''Convert seq (string with DNA alphabet) to numpy array with integers 0,1,2,3.'''
    return np.array([_base2bits[base] for base in seq],dtype=np.ubyte)

def seqs2arr(seqs):
    '''Return numpy 2D array converting the given DNA sequences to integers.'''
    seqLen = len(seqs[0])
    for seq in seqs:
        if len(seq) != seqLen:
            raise ValueError('All sequences in seqs must be equal length')
    numSeqs = len(seqs)
    arr = np.empty((numSeqs, seqLen), dtype=np.ubyte)
    for i in range(numSeqs):
        arr[i] = [_base2bits[base] for base in seqs[i]]
    return arr

def arr2seq(arr):
    basesCh = map(lambda base: _bits2base[base], arr)
    return ''.join(basesCh)

def makeArrayWithAllDNASeqs(length,bases=('A','C','G','T')):
    '''Return 2D numpy array with all DNA sequences of given length in
    lexicographic order. Bases contains bases to be used: ('A','C','G','T') by
    default, but can be set to a subset of these.

    Uses the encoding described in the documentation for DNASeqList. The result
    is a 2D array, where each row represents a DNA sequence, and that row
    has one byte per base.'''

    if not set(bases) <= set(('A','C','G','T')):
        raise ValueError("bases must be a subset of ('A','C','G','T'); cannot be %s" % bases)
    if len(bases) == 0:
        raise ValueError('bases cannot be empty')

    numbases = len(bases)
    numseqs = numbases**length

#     shift = np.arange(2*(length-1), -1, -2)
#     nums = np.repeat(np.arange(numseqs), length)
#     nums2D = nums.reshape([numseqs, length])
#     shifts = np.tile(0b11 << shift, numseqs)
#     shifts2D = shifts.reshape([numseqs, length])
#     arr = (shifts2D & nums2D) >> shift

    # the former code took up too much memory (using int32 or int64)
    # the following code makes sure it's 1 byte per base
    powers_numbases = [numbases**k for k in range(length)]
#         bases = np.array([_base2bits['A'], _base2bits['C'], _base2bits['G'], _base2bits['T']], dtype=np.ubyte)
    basebits = map(lambda base: _base2bits[base], bases)
    bases = np.array(basebits, dtype=np.ubyte)

    list_of_arrays = False
    if list_of_arrays:
        # this one seems to be faster but takes more memory, probably because
        # just before the last command there are two copies of the array
        # in memory at once
        columns = []
        for i,j,c in zip(reversed(powers_numbases), powers_numbases, range(length)):
            columns.append(np.tile(np.repeat(bases, i), j))
        arr = np.vstack(columns).transpose()
    else:
        # this seems to be slightly slower but takes less memory, since it
        # allocates only the final array, plus one extra column of that
        # array at a time
        arr = np.empty((numseqs, length), dtype=np.ubyte)
        for i,j,c in zip(reversed(powers_numbases), powers_numbases, range(length)):
            arr[:,c] = np.tile(np.repeat(bases, i), j)

    return arr

# from http://www.bogotobogo.com/python/python_longest_common_substring_lcs_algorithm_generalized_suffix_tree.php
@lru_cache(maxsize=10000000)
def lcs(S,T):
    '''Longest common substring of two Python strings.'''
    m = len(S)
    n = len(T)
    counter = [[0]*(n+1) for x in range(m+1)]
    longest = 0
    lcs_set = set()
    for i in range(m):
        for j in range(n):
            if S[i] == T[j]:
                c = counter[i][j] + 1
                counter[i+1][j+1] = c
                if c > longest:
                    lcs_set = set()
                    longest = c
                    lcs_set.add(S[i-c+1:i+1])
                elif c == longest:
                    lcs_set.add(S[i-c+1:i+1])

    return lcs_set

#@lru_cache(maxsize=10000000)
def longest_common_substring(a1, a2, vectorized=True):
    '''Return start and end indices (a1start, a2start, length) of longest common
    substring (subarray) of 1D arrays a1 and a2.'''
    assert len(a1.shape) == 1
    assert len(a2.shape) == 1
    counter = np.zeros(shape=(len(a1)+1,len(a2)+1), dtype=np.int)
    a1idx_longest = a2idx_longest = -1
    len_longest = 0

    if vectorized:
        for i1 in range(len(a1)):
            idx = (a2 == a1[i1])
            idx_shifted = np.hstack([[False], idx])
            counter[i1+1,idx_shifted] = counter[i1,idx]+1
        idx_longest = np.unravel_index(np.argmax(counter), counter.shape)
        if idx_longest[0] > 0:
            len_longest = counter[idx_longest]
            a1idx_longest = int(idx_longest[0] - len_longest)
            a2idx_longest = int(idx_longest[1] - len_longest)
    else:
        for i1 in range(len(a1)):
            for i2 in range(len(a2)):
                if a1[i1] == a2[i2]:
                    c = counter[i1,i2] + 1
                    counter[i1+1,i2+1] = c
                    if c > len_longest:
                        len_longest = c
                        a1idx_longest = i1+1-c
                        a2idx_longest = i2+1-c
    return (a1idx_longest, a2idx_longest, len_longest)

#@lru_cache(maxsize=10000000)
def longest_common_substrings_singlea1(a1, a2s):
    '''Return start and end indices (a1starts, a2starts, lengths) of longest common
    substring (subarray) of 1D array a1 and rows of 2D array a2s.

    If length[i]=0, then a1starts[i]=a2starts[i]=0 (not -1), so be sure to check
    length[i] to see if any substrings actually matched.'''
    assert len(a1.shape) == 1
    assert len(a2s.shape) == 2
    numa2s = a2s.shape[0]
    len_a1 = len(a1)
    len_a2 = a2s.shape[1]
    counter = np.zeros(shape=(len_a1+1, numa2s, len_a2+1), dtype=np.int)

    for i1 in range(len(a1)):
        idx = (a2s == a1[i1])
        idx_shifted = np.insert(idx, 0, np.zeros(numa2s, dtype=np.bool), axis=1)
        counter[i1+1,idx_shifted] = counter[i1,idx]+1

    counter = np.swapaxes(counter, 0, 1)

    counter_flat = counter.reshape(numa2s, (len_a1+1)*(len_a2+1))
    idx_longest_raveled = np.argmax(counter_flat, axis=1)
    len_longest = counter_flat[np.arange(counter_flat.shape[0]), idx_longest_raveled]

    idx_longest = np.unravel_index(idx_longest_raveled, dims=(len_a1+1, len_a2+1))
    a1idx_longest = idx_longest[0] - len_longest
    a2idx_longest = idx_longest[1] - len_longest

    return (a1idx_longest, a2idx_longest, len_longest)

#@lru_cache(maxsize=10000000)
def longest_common_substrings_product(a1s, a2s):
    '''Return start and end indices (a1starts, a2starts, lengths) of longest common
    substring (subarray) of each pair in the cross product of rows of a1s and a2s.

    If length[i]=0, then a1starts[i]=a2starts[i]=0 (not -1), so be sure to check
    length[i] to see if any substrings actually matched.'''
    numa1s = a1s.shape[0]
    numa2s = a2s.shape[0]
    a1s_cp = np.repeat(a1s, numa2s, axis=0)
    a2s_cp = np.tile(a2s, (numa1s, 1))

    a1idx_longest, a2idx_longest, len_longest = _longest_common_substrings_pairs(a1s_cp, a2s_cp)

    a1idx_longest = a1idx_longest.reshape(numa1s, numa2s)
    a2idx_longest = a2idx_longest.reshape(numa1s, numa2s)
    len_longest = len_longest.reshape(numa1s, numa2s)

    return (a1idx_longest, a2idx_longest, len_longest)

def pair_index(n):
    index = np.fromiter(it.chain.from_iterable(it.combinations(xrange(n), 2)), int, count=n*(n-1))
    return index.reshape(-1, 2)

def _longest_common_substrings_pairs(a1s, a2s):
    assert len(a1s.shape) == 2
    assert len(a2s.shape) == 2
    assert a1s.shape[0] == a2s.shape[0]

    numpairs = a1s.shape[0]

    len_a1 = a1s.shape[1]
    len_a2 = a2s.shape[1]

    counter = np.zeros(shape=(len_a1+1, numpairs, len_a2+1), dtype=np.int)

    for i1 in range(len_a1):
        a1s_cp_col = a1s[:,i1].reshape(numpairs,1)
        a1s_cp_col_rp = np.repeat(a1s_cp_col, len_a2, axis=1)

        idx = (a2s == a1s_cp_col_rp)
        idx_shifted = np.hstack([np.zeros(shape=(numpairs,1), dtype=np.bool), idx])
        counter[i1+1,idx_shifted] = counter[i1,idx]+1

    counter = np.swapaxes(counter, 0, 1)

    counter_flat = counter.reshape(numpairs, (len_a1+1)*(len_a2+1))
    idx_longest_raveled = np.argmax(counter_flat, axis=1)
    len_longest = counter_flat[np.arange(counter_flat.shape[0]), idx_longest_raveled]

    idx_longest = np.unravel_index(idx_longest_raveled, dims=(len_a1+1, len_a2+1))
    a1idx_longest = idx_longest[0] - len_longest
    a2idx_longest = idx_longest[1] - len_longest

    return (a1idx_longest, a2idx_longest, len_longest)

def longest_common_substrings_all_pairs_strings(seqs1, seqs2):
    '''For Python strings'''
    a1s = seqs2arr(seqs1)
    a2s = seqs2arr(seqs2)
    return _longest_common_substrings_pairs(a1s, a2s)


def _strongest_common_substrings_all_pairs_return_energies_and_counter(a1s, a2s, T):
    assert len(a1s.shape) == 2
    assert len(a2s.shape) == 2
    assert a1s.shape[0] == a2s.shape[0]

    numpairs = a1s.shape[0]
    len_a1 = a1s.shape[1]
    len_a2 = a2s.shape[1]
    counter = np.zeros(shape=(len_a1+1, numpairs, len_a2+1), dtype=np.int)
    energies = np.zeros(shape=(len_a1+1, numpairs, len_a2+1), dtype=np.float)

#     if not loop_energies:
    loop_energies = calculate_loop_energies(T)

    prev_match_idxs = prev_match_shifted_idxs = None

    for i1 in range(len_a1):
        a1s_col = a1s[:,i1].reshape(numpairs,1)
        a1s_col_rp = np.repeat(a1s_col, len_a2, axis=1)

        # find matching chars and extend length of substring
        match_idxs = (a2s == a1s_col_rp)
        match_shifted_idxs = np.hstack([np.zeros(shape=(numpairs,1), dtype=np.bool), match_idxs])
        counter[i1+1,match_shifted_idxs] = counter[i1,match_idxs] + 1

        if i1 > 0:
            # calculate energy if matching substring has length > 1
            prev_bases = a1s[:,i1-1]
            cur_bases = a1s[:,i1]
            loops = (prev_bases << 2) + cur_bases
            latest_energies = loop_energies[loops].reshape(numpairs, 1)
            latest_energies_rp = np.repeat(latest_energies, len_a2, axis=1)
            match_idxs_false_at_end = np.hstack([match_idxs, np.zeros(shape=(numpairs,1), dtype=np.bool)])
            both_match_idxs = match_idxs_false_at_end & prev_match_shifted_idxs
            prev_match_shifted_shifted_idxs = np.hstack([np.zeros(shape=(numpairs,1), dtype=np.bool), prev_match_shifted_idxs])[:,:-1]
            both_match_shifted_idxs = match_shifted_idxs & prev_match_shifted_shifted_idxs
            energies[i1+1,both_match_shifted_idxs] = energies[i1,both_match_idxs] + latest_energies_rp[both_match_idxs]

#         prev_match_idxs = match_idxs
        prev_match_shifted_idxs = match_shifted_idxs

    counter = counter.swapaxes(0, 1)
    energies = energies.swapaxes(0, 1)

    return counter, energies

def weighted_energies_common_substrings(seqs1, seqs2, T):
    a1s = seqs2arr(seqs1)
    a2s = seqs2arr(seqs2)
    numpairs = a1s.shape[0]
    len_a1 = a1s.shape[1]
    len_a2 = a2s.shape[1]

    counter, energies = _strongest_common_substrings_all_pairs_return_energies_and_counter(a1s, a2s, T)

    counter_shift_upleft = np.roll(np.roll(counter, -1, axis=1), -1, axis=2)
    not_maximal_substring_idxs = (counter == 0) | (counter_shift_upleft != 0)
    energies[not_maximal_substring_idxs] = 0
    energies_flat = energies.reshape((numpairs, (len_a1+1)*(len_a2+1)))
    weighted_energies = np.zeros(numpairs, dtype=np.float)

    raise NotImplementedError()

    return list(weighted_energies)

import math

def internal_loop_penalty(n, T):
    return 1.5 + (2.5*0.002*T*math.log(1+n))

def _mfes_array(a1s, a2s, T):
    assert len(a1s.shape) == 2
    assert len(a2s.shape) == 2
    assert a1s.shape[0] == a2s.shape[0]

    numpairs = a1s.shape[0]
    len_a1 = a1s.shape[1]
    len_a2 = a2s.shape[1]
    d = np.zeros(shape=(len_a1+1, len_a2+1, numpairs), dtype=np.int)
    energies = np.zeros(shape=(len_a1+1, len_a2+1, numpairs), dtype=np.float)

#     if not loop_energies:
    loop_energies = calculate_loop_energies(T)

    for i1 in range(len_a1):
        for i2 in range(len_a2):
            recursive_energy_1 = energies[i1-1,i2,:] - internal_loop_penalty( d[i1-1,i2,:] + 1, T )
            recursive_energy_2 = energies[i1,i2-1,:] - internal_loop_penalty( d[i1,i2-1,:] + 1, T )
            match_idxs = (a1s[i1,i2,:] == a2s[i1,i2,:])


    return energies

def mfes(seqs1, seqs2, T):
    a1s = seqs2arr(seqs1)
    a2s = seqs2arr(seqs2)
    a2swc = wc(a2s)

    energies = _mfes_array(a1s, a2swc, T)

    return list(energies)

def _strongest_common_substrings_all_pairs(a1s, a2s, T):
    numpairs = a1s.shape[0]
    len_a1 = a1s.shape[1]
    len_a2 = a2s.shape[1]

    counter, energies = _strongest_common_substrings_all_pairs_return_energies_and_counter(a1s, a2s, T)

    counter_flat = counter.reshape(numpairs, (len_a1+1)*(len_a2+1))
    energies_flat = energies.reshape(numpairs, (len_a1+1)*(len_a2+1))

    idx_strongest_raveled = np.argmax(energies_flat, axis=1)
    len_strongest = counter_flat[np.arange(counter_flat.shape[0]), idx_strongest_raveled]
    energy_strongest = energies_flat[np.arange(counter_flat.shape[0]), idx_strongest_raveled]

    idx_strongest = np.unravel_index(idx_strongest_raveled, dims=(len_a1+1, len_a2+1))
    a1idx_strongest = idx_strongest[0] - len_strongest
    a2idx_strongest = idx_strongest[1] - len_strongest

    return (a1idx_strongest, a2idx_strongest, len_strongest, energy_strongest)


def strongest_common_substrings_all_pairs_string(seqs1, seqs2, T):
    '''For Python strings representing DNA; checks for reverse complement matches
    rather than direct matches, and evaluates nearest neighbor energy, returning
    indices lengths, and energies of strongest complementary substrings.'''
    a1s = seqs2arr(seqs1)
    a2s = seqs2arr(seqs2)
    a1idx_strongest, a2idx_strongest, len_strongest, energy_strongest = _strongest_common_substrings_all_pairs(a1s, wc(a2s), T)
    return (list(a1idx_strongest), list(a2idx_strongest), list(len_strongest), list(energy_strongest))

def energies_strongest_common_substrings(seqs1, seqs2, T):
    a1s = seqs2arr(seqs1)
    a2s = seqs2arr(seqs2)
    a1idx_strongest, a2idx_strongest, len_strongest, energy_strongest = _strongest_common_substrings_all_pairs(a1s, wc(a2s), T)
    return list(energy_strongest)


def binding_heuristic(seq1, seq2, num_subs=3, len_subs=4):
    """Returns True if there are num_subs, each of length len_subs, complementary
    between seq1 and seq2, which are Python strings."""
    all_subs = makeArrayWithAllDNASeqs(len_subs)
    raise NotImplementedError()

def binding_heuristic_numpy(seq, seqs, num_subs=3, len_subs=4):
    """Returns array of booleans indicating whether
    binding_heuristic(seq,seq[i],num_subs,len_subs) for each seq[i] in seqs (a
    DNASeqList), where seq is a numpy array of integers from [0,1,2,3] or
    a Python string from ['A','C','G','T']"""
    raise NotImplementedError()

class DNASeqList(object):
    '''Represents a list of DNA sequences of identical length.

    Uses a (noncompact) internal representation using 8 bits per base, stored
    in a numpy 2D array of bytes.

    The code used is A-->0, C-->1, G-->2, T-->3'''

    def __init__(self,length=0,alphabet=('A','C','G','T'),seqs=[],seqarr=None,
                 fileName=None,binaryFileName=None):
        '''Creates a set of DNA sequences.

        Create either all sequences of a given length if seqs is not specified,
        or all sequences in seqs if seqs is specified. If neither is specified
        then all sequences of length 3 are created.'''
        if seqarr is not None:
            self.seqarr = seqarr
            self.numseqs, self.seqlen = seqarr.shape
        elif seqs:
            self.seqlen = len(seqs[0])
            for seq in seqs:
                if len(seq) != self.seqlen:
                    raise ValueError('All sequences in seqs must be equal length')
            self.numseqs = len(seqs)
            self.seqarr = seqs2arr(seqs)
        elif fileName is not None:
            self._readFromFile(fileName)
        elif binaryFileName is not None:
            self._readFromBinaryFile(binaryFileName)
        elif length:
            self.seqlen = length
            self.seqarr = makeArrayWithAllDNASeqs(self.seqlen, alphabet)
            self.numseqs = len(self.seqarr)
        else:
            raise ValueError('at least one of length, seqs, seqarr, fileName, or binaryFileName must be specified')
        self.shift = np.arange(2*(self.seqlen-1), -1, -2)

    def __len__(self):
        return self.numseqs

    def __contains__(self, seq):
        if len(seq) != self.seqlen:
            return False
        arr = seq2arr(seq)
        return np.any(~np.any(self.seqarr - arr, 1))

    def _readFromFile(self,fileName):
        '''Reads from fileName in the format defined in writeToFile.
        Only meant to be called from constructor.'''
        with open(fileName,'r+') as f:
            firstLine = f.readline()
            numSeqsStr,seqLenStr,temperature = firstLine.split()
            self.numseqs = int(numSeqsStr)
            self.seqlen = int(seqLenStr)
            self.seqarr = np.empty((self.numseqs, self.seqlen), dtype=np.ubyte)
            for i in range(self.numseqs):
                line = f.readline()
                seq = line.strip()
                self.seqarr[i] = [_base2bits[base] for base in seq]

    def writeToFile(self,fileName):
        '''Writes text file describing DNA sequence list, in format

        numseqs seqlen
        seq1
        seq2
        seq3
        ...

        where numseqs, seqlen are integers, and seq1,
        ... are strings from {A,C,G,T}'''
        with open(fileName,'w+') as f:
            f.write(str(self.numseqs) + ' ' + str(self.seqlen) + '\n')
            for i in range(self.numseqs):
                f.write(self.getSeqStr(i) + '\n')

    def _readFromBinaryFile(self,fileName):
        '''Reads from fileName in the format defined in writeToBinaryFile.
        Only meant to be called from constructor.'''
        f = np.load(fileName)
        self.numseqs = f['numseqs']
        self.seqlen = f['seqlen']
        packedbits = f['packedbits']
        self.seqarr = unpack(packedbits,self.seqlen)

    def writeToBinaryFile(self,fileName):
        '''Writes binary file describing DNA sequence list, in format used by
        numpy.savez to store multiple arrays in a zip file, using a packed bits
        representation for the DNA.'''
        packedbits = pack(self.seqarr)
        np.savez(fileName, numseqs=self.numseqs,
                           seqlen=self.seqlen,
                           packedbits=packedbits)

    def _readFromCompressedBinaryFile(self,fileName):
        '''Reads from fileName in the format defined in writeToCompressedBinaryFile.
        Only meant to be called from constructor.'''
        f = np.load(fileName)
        self.numseqs = f['numseqs']
        self.seqlen = f['seqlen']
        self.seqarr = f['seqarr']

    def writeToCompressedBinaryFile(self,fileName):
        '''Writes binary file describing DNA sequence list, in format used by
        numpy.savez_compressed to store multiple arrays in a zip file.'''
        np.savez_compressed(fileName, numseqs=self.numseqs,
                                      seqlen=self.seqlen,
                                      seqarr=self.seqarr)

    def wcenergy(self,idx,temperature):
        '''Return energy of idx'th sequence binding to its complement.'''
        return wcenergy(self.seqarr[idx],temperature)

    def __repr__(self):
        return 'DNASeqSet(seqs={})'.format(str([self[i] for i in range(self.numseqs)]))

    def __str__(self):
        if self.numseqs <= 64:
            ret = [self.getSeqStr(i) for i in range(self.numseqs)]
            return ','.join(ret)
        else:
            ret = [self.getSeqStr(i) for i in range(3)] + ['...'] + \
                  [self.getSeqStr(i) for i in range(self.numseqs-3,self.numseqs)]
            return ','.join(ret)

    def shuffle(self):
        import numpy.random as nprand
        nprand.shuffle(self.seqarr)

    def toList(self):
        '''Return list of strings representing the sequences, e.g. ['ACG','TAA']'''
        return [self.getSeqStr(idx) for idx in range(self.numseqs)]

    def getSeqStr(self,idx):
        '''Return idx'th DNA sequence as a string.'''
        return arr2seq(self.seqarr[idx])

    def getSeqsStrList(self,slice):
        '''Return a list of strings specified by slice.'''
        bases_lst = self.seqarr[slice]
        ret = []
        for bases in bases_lst:
            basesCh = map(lambda base: _bits2base[base], bases)
            ret.append(''.join(basesCh))
        return ret

    def __getitem__(self,idx):
        if type(idx) is int:
            return self.getSeqStr(idx)
        elif type(idx) is slice:
            return self.getSeqsStrList(idx)
        else:
            raise ValueError('idx must be int or slice')

    def pop(self):
        '''Remove and return last seq, as a string.'''
        seqStr = self.getSeqStr(-1)
        self.seqarr = np.delete(self.seqarr, -1, 0)
        self.numseqs -= 1
        return seqStr

    def pop_array(self):
        '''Remove and return last seq, as a string.'''
        arr = self.seqarr[-1]
        self.seqarr = np.delete(self.seqarr, -1, 0)
        self.numseqs -= 1
        return arr

    def append_seq(self, newseq):
        self.append_arr(seq2arr(newseq))

    def append_arr(self, newarr):
        self.seqarr = np.vstack([self.seqarr, newarr])
        self.numseqs += 1

    def filter_hamming(self, threshold):
        seq = self.pop_array()
        arr_keep = np.array([seq])
        self.shuffle()
        while self.seqarr.shape[0] > 0:
            seq = self.pop_array()
            while self.seqarr.shape[0] > 0:
                hamming_min = np.min(np.sum(np.bitwise_xor(arr_keep, seq) != 0, axis=1))
                too_close = (hamming_min < threshold)
                if not too_close:
                    break
                seq = self.pop_array()
            arr_keep = np.vstack([arr_keep, seq])
        self.seqarr = arr_keep
        self.numseqs = self.seqarr.shape[0]

    def hamming_min(self, arr):
        '''Returns minimum Hamming distance between arr and any sequence in
        this DNASeqList.'''
        distances = np.sum(np.bitwise_xor(self.seqarr, arr) != 0, axis=1)
        return np.min(distances)


#     def getSeqLong(self,idx):
#         '''Return idx'th DNA sequence as a long in bit-compressed format.'''
#         seq = self.seqarr[idx]
#         shiftedDigits = seq << self.shift
#         lng = np.sum(shiftedDigits)
#         return lng

    def filter_energy(self, low, high, temperature):
        '''Return new DNASeqList with seqs whose wc complement energy is within
        the given range.'''
        wcenergies = calculateWCenergies(self.seqarr, temperature)
        within_range = (low <= wcenergies) & (wcenergies <= high)
        new_seqarr = self.seqarr[within_range]
        return DNASeqList(seqarr=new_seqarr)

    def energies(self, temperature):
        wcenergies = calculateWCenergies(self.seqarr, temperature)
        return wcenergies

    def filter_PF(self, low, high, temperature):
        '''Return new DNASeqList with seqs whose wc complement energy is within
        the given range, according to NUPACK.'''
        raise NotImplementedError('this was assuming energies get stored, which no longer is the case')
        pfenergies = np.zeros(self.wcenergies.shape)
        i = 0
        print 'searching through %d sequences for PF energies' % self.numseqs
        for seq in self.toList():
            energy = sst_dsd.duplex(seq, temperature)
            pfenergies[i] = energy
            i += 1
            if i % 100 == 0:
                print 'searched %d so far' % i
        within_range = (low <= pfenergies) & (pfenergies <= high)
        new_seqarr = self.seqarr[within_range]
        return DNASeqList(seqarr=new_seqarr)

    def filter_endGC(self):
        '''Remove any sequence with A or T on the end. Also remove domains that
        do not have an A or T either next to that base, or one away. Otherwise
        we could get a domain ending in {C,G}^3, which, placed next to any
        domain ending in C or G, will create a substring in {C,G}^4 and be
        rejected if we are filtering those.'''
        left = self.seqarr[:,0]
        right = self.seqarr[:,-1]
        leftP1 = self.seqarr[:,1]
        leftP2 = self.seqarr[:,2]
        rightM1 = self.seqarr[:,-2]
        rightM2 = self.seqarr[:,-3]
        cbits = _base2bits['C']
        gbits = _base2bits['G']
        abits = _base2bits['A']
        tbits = _base2bits['T']
        good = ( ((left == cbits) | (left == gbits)) & ((right == cbits) | (right == gbits)) &
                 ((leftP1 == abits) | (leftP1 == tbits) | (leftP2 == abits) | (leftP2 == tbits)) &
                 ((rightM1 == abits) | (rightM1 == tbits) | (rightM2 == abits) | (rightM2 == tbits)) )
        seqarrpass = self.seqarr[good]
        return DNASeqList(seqarr=seqarrpass)

    def filter_endAT(self, gc_near_end=False):
        '''Remove any sequence with C or G on the end. Also, if gc_near_end is True,
        remove domains that do not have an C or G either next to that base,
        or one away, to prevent breathing.'''
        left = self.seqarr[:,0]
        right = self.seqarr[:,-1]
        abits = _base2bits['A']
        tbits = _base2bits['T']
        good = ((left == abits) | (left == tbits)) & ((right == abits) | (right == tbits))
        if gc_near_end:
            cbits = _base2bits['C']
            gbits = _base2bits['G']
            leftP1 = self.seqarr[:,1]
            leftP2 = self.seqarr[:,2]
            rightM1 = self.seqarr[:,-2]
            rightM2 = self.seqarr[:,-3]
            good = ( good &
                   ((leftP1 == cbits) | (leftP1 == gbits) | (leftP2 == cbits) | (leftP2 == gbits)) &
                   ((rightM1 == cbits) | (rightM1 == gbits) | (rightM2 == cbits) | (rightM2 == gbits)) )
        seqarrpass = self.seqarr[good]
        return DNASeqList(seqarr=seqarrpass)

    def filter_base_nowhere(self, base):
        '''Remove any sequence that has given base anywhere.'''
        good = (self.seqarr != _base2bits[base]).all(axis=1)
        seqarrpass = self.seqarr[good]
        return DNASeqList(seqarr=seqarrpass)

    def filter_base_count(self, base, low, high):
        '''Remove any sequence not satisfying low <= #base <= high.'''
        sumarr = np.sum(self.seqarr == _base2bits[base], axis=1)
        good = (low <= sumarr) & (sumarr <= high)
        seqarrpass = self.seqarr[good]
        return DNASeqList(seqarr=seqarrpass)

    def filter_base_at_pos(self, pos, base):
        '''Remove any sequence that does not have given base at position pos.'''
        mid = self.seqarr[:,pos]
        good = (mid == _base2bits[base])
        seqarrpass = self.seqarr[good]
        return DNASeqList(seqarr=seqarrpass)

    def filter_substring(self,subs):
        '''Remove any sequence with any elements from subs as a substring.'''
        if len(set([len(sub) for sub in subs])) != 1:
            raise ValueError('All substrings in subs must be equal length: %s' % subs)
        sublen=len(subs[0])
        subints = [[_base2bits[base] for base in sub] for sub in subs]
        powarr=[4**k for k in range(sublen)]
        subvals = np.dot(subints,powarr)
        toeplitz = create_toeplitz(self.seqlen,len(sub))
        convolution = np.dot(toeplitz, self.seqarr.transpose())
        passall = np.ones(self.numseqs,dtype=np.bool)
        for subval in subvals:
            passsub = np.all(convolution != subval,axis=0)
            passall = passall & passsub
        seqarrpass = self.seqarr[passall]
        return DNASeqList(seqarr=seqarrpass)

    def filterSeqsByGQuad(self,fast=True):
        '''Removes any sticky ends with 4 G's in a row (a G-quadruplex).'''
        return self.filter_substring(['GGGG'])

    def filterSeqsByGQuadCQuad(self,fast=True):
        '''Removes any sticky ends with 4 G's or C's in a row (a quadruplex).'''
        return self.filter_substring(['GGGG','CCCC'])

def create_toeplitz(seqlen,sublen):
    '''Creates a toeplitz matrix, useful for finding subsequences.'''
    powarr = [4**k for k in range(sublen)]
    if _SCIPY_AVAIL:
        import scipy.linalg as linalg
        toeplitz = linalg.toeplitz([1]+[0]*(seqlen-sublen),
                                   powarr+[0]*(seqlen-sublen))
    else: #This is a fix for the fact that numpypy doesn't have scipy.
        rows = seqlen-(sublen-1)
        cols = seqlen
        toeplitz = np.zeros((rows,cols),dtype=np.int)
        toeplitz[:,0:sublen] = [powarr]*rows
        shift = range(rows)
        for i in range(rows):
            toeplitz[i] = np.roll(toeplitz[i],shift[i])
    return toeplitz

from lru_cache import lru_cache
@lru_cache(maxsize=32)
def calculate_loop_energies(temperature):
    '''Get SantaLucia and Hicks nearest-neighbor loop energies for given temperature,
    1 M Na+. '''
    return -(_dH - (temperature+273.15)*_dS/1000.0)
    # SantaLucia & Hicks' values are in cal/mol/K for dS, and kcal/mol for dH.
    # Here we divide dS by 1000 to get the RHS term into units of kcal/mol/K
    # which gives an overall dG in units of kcal/mol.
    # One reason we want dG to be in units of kcal/mol is to
    # give reasonable/readable numbers close to 0 for dG(Assembly).
    # The reason we flip the sign is that, by convention, in the kTAM, G_se
    # (which is computed from the usually negative dG here) is usually positive.


# _dH and _dS come from Table 1 in SantaLucia and Hicks, Annu Rev Biophys Biomol Struct. 2004;33:415-40.
#                 AA   AC   AG   AT   CA   CC    CG   CT
_dH = np.array([-7.6,-8.4,-7.8,-7.2,-8.5,-8.0,-10.6,-7.8,
                # GA   GC   GG   GT   TA   TC   TG   TT
                -8.2,-9.8,-8.0,-8.4,-7.2,-8.2,-8.5,-7.6],
               dtype=np.float32)

#                  AA    AC    AG    AT    CA    CC    CG    CT
_dS = np.array([-21.3,-22.4,-21.0,-20.4,-22.7,-19.9,-27.2,-21.0,
                #  GA    GC    GG    GT    TA    TC    TG    TT
                -22.2,-24.4,-19.9,-22.4,-21.3,-22.2,-22.7,-21.3],
               dtype=np.float32)


#  AA  AC  AG  AT  CA  CC  CG  CT  GA  GC  GG  GT  TA  TC  TG  TT
#  00  01  02  03  10  11  12  13  20  21  22  23  30  31  32  34

# nearest-neighbor energies for Watson-Crick complements at 37C
# (Table 1 in SantaLucia and Hicks 2004)
# ordering of array is
# #                      AA    AC    AG    AT    CA    CC    CG    CT
# _nndGwc = np.array([-1.00,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,
#                     #  GA    GC    GG    GT    TA    TC    TG    TT
#                     -1.30,-2.24,-1.84,-1.44,-0.58,-1.30,-1.45,-1.00],
#                    dtype=np.float32)
#                    # AA   AC   AG   AT   CA   CC   CG   CT
#_nndGwc = np.array([1.00,1.44,1.28,0.88,1.45,1.84,2.17,1.28,
#                    # GA   GC   GG   GT   TA   TC   TG   TT
#                    1.30,2.24,1.84,1.44,0.58,1.30,1.45,1.00],
#                   dtype=np.float32)
#_nndGwcStr = {'AA':1.00,'AC':1.44,'AG':1.28,'AT':0.88,'CA':1.45,'CC':1.84,
#              'CG':2.17,'CT':1.28,'GA':1.30,'GC':2.24,'GG':1.84,'GT':1.44,
#              'TA':0.58,'TC':1.30,'TG':1.45,'TT':1.00}

# nearest-neighbor energies for single mismatches (Table 2 in SantaLucia)
# ordering of array is
#                  # GA/CA GA/CG    AG    AT    CA    CC    CG    CT   GA    GC    GG    GT    TA   TC    TG   TT
#_nndGsmm = np.array([0.17,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.3,-2.24,-1.84,-1.44,-0.58,-1.3,-1.45,-1.0], dtype=np.float32)

_all_pairs = [((i<<2)+j,_bits2base[i]+_bits2base[j])
              for i in range(4) for j in range(4)]

@lru_cache(maxsize=32)
def calculate_loop_energies_dict(temperature):
    loop_energies = calculate_loop_energies(temperature)
    return {pair[1]: loop_energies[pair[0]] for pair in _all_pairs}

@lru_cache(maxsize=100000)
def wcenergy(seq,temperature):
    '''Return the wc energy of seq binding to its complement.'''
    loop_energies = calculate_loop_energies_dict(temperature)
    return sum(loop_energies[seq[i:i+2]] for i in range(len(seq)-1))

def wcenergiesStr(seqs,temperature):
    seqarr = seqs2arr(seqs)
    return list(calculateWCenergies(seqarr,temperature))

def wcenergyStr(seq,temperature):
    seqarr = seqs2arr([seq])
    return list(calculateWCenergies(seqarr,temperature))[0]

def hash_ndarray(arr):
    writeable = arr.flags.writeable
    if writeable:
        arr.flags.writeable = False
    h = hash(arr.data)
    arr.flags.writeable = writeable
    return h

CACHE_WC = True
_calculateWCenergies_cache = None
_calculateWCenergies_cache_hash = 0

def calculateWCenergies(seqarr,temperature):
    '''Calculate and store in an array all energies of all sequences in seqarr
    with their Watson-Crick complements.'''
    global _calculateWCenergies_cache
    global _calculateWCenergies_cache_hash
    if CACHE_WC and _calculateWCenergies_cache is not None:
        if _calculateWCenergies_cache_hash == hash_ndarray(seqarr):
            return _calculateWCenergies_cache
    loop_energies = calculate_loop_energies(temperature)
    leftIndexBits = seqarr[:,:-1] << 2
    rightIndexBits = seqarr[:,1:]
    pairIndices = leftIndexBits + rightIndexBits
    pairEnergies = loop_energies[pairIndices]
    if CACHE_WC:
        _calculateWCenergies_cache = np.sum(pairEnergies, axis=1)
        _calculateWCenergies_cache_hash = hash_ndarray(_calculateWCenergies_cache)
    return _calculateWCenergies_cache

def wc(seqarr):
    '''Return array of complements of sequences in seqarr'''
    return (3 - seqarr)[:,::-1]














def filterByCommonDensityLinear(seqs1, seqs2, delta=0.01, numPoints=100):
    '''Find a point where seqs1 and seqs2 both have high density and filter
    everything not within delta fraction of that point.

    Assumes both sequence lists have been sorted by binding energy.'''
    minPoint = np.min([seqs1.wcenergies[0], seqs2.wcenergies[0]])
    maxPoint = np.max([seqs1.wcenergies[-1], seqs2.wcenergies[-1]])
    points = np.linspace(minPoint, maxPoint, numPoints, endpoint=False)
    posDensest = 0
    numCloseMaxMin = 0
    for i in range(numPoints):
        low = points[i]
        diff = delta*low
        high = low+diff
        posLow1  = np.searchsorted(seqs1.wcenergies, low, side='left')
        posHigh1 = np.searchsorted(seqs1.wcenergies, high, side='right')
        posLow2  = np.searchsorted(seqs2.wcenergies, low, side='left')
        posHigh2 = np.searchsorted(seqs2.wcenergies, high, side='right')
        numClose1 = posHigh1 - posLow1
        numClose2 = posHigh2 - posLow2
        numCloseMin = np.min([numClose1, numClose2])
        if numCloseMaxMin < numCloseMin:
            numCloseMaxMin = numCloseMin
            posDensest = i
    low = points[posDensest]
    diff = delta*low
    high = low + diff
    filteredSeqs1 = seqs1.filter_energy(low, high)
    filteredSeqs2 = seqs2.filter_energy(low, high)
    return (filteredSeqs1, filteredSeqs2)

def filterByCommonDensityBinary(seqs1, seqs2, delta=0.01):
    '''Find a point where seqs1 and seqs2 both have high density and filter
    everything not within delta fraction of that point.

    Conducts a binary search, starting with low and high at the endpoints of
    the energy ranges spanned by both lists, and stops when
    (high-low) <= delta*low (i.e., the largest difference in energies is
    less than delta fraction of low)

    Assumes both sequence lists have been sorted by binding energy.'''
    low = np.min([seqs1.wcenergies[0], seqs2.wcenergies[0]])
    high = np.max([seqs1.wcenergies[-1], seqs2.wcenergies[-1]])
    mid = (low + high) / 2
    posLow1  = np.searchsorted(seqs1.wcenergies, low, side='left')
    posLow2  = np.searchsorted(seqs2.wcenergies, low, side='left')
    posHigh1 = np.searchsorted(seqs1.wcenergies, high, side='right')
    posHigh2 = np.searchsorted(seqs2.wcenergies, high, side='right')
    while high - low > delta*low:
        mid = (low + high) / 2
        posMid1 = np.searchsorted(seqs1.wcenergies, mid, side='right')
        posMid2 = np.searchsorted(seqs2.wcenergies, mid, side='right')
        leftCount1 = posMid1 - posLow1
        rightCount1 = posHigh1 - posMid1
        leftCount2 = posMid2 - posLow2
        rightCount2 = posHigh2 - posMid2
        leftMinCount = np.min([leftCount1, leftCount2])
        rightMinCount = np.min([rightCount1, rightCount2])
        if leftMinCount > rightMinCount: # go left
            high = mid
            posHigh1 = np.searchsorted(seqs1.wcenergies, high, side='right')
            posHigh2 = np.searchsorted(seqs2.wcenergies, high, side='right')
        else: # go right
            low = mid
            posLow1  = np.searchsorted(seqs1.wcenergies, low, side='left')
            posLow2  = np.searchsorted(seqs2.wcenergies, low, side='left')
    high = low * (1 + delta)
    filteredSeqs1 = seqs1.filter_energy(low, high)
    filteredSeqs2 = seqs2.filter_energy(low, high)
    return (filteredSeqs1, filteredSeqs2)

def filterByCommonDensityVectorized(seqs1, seqs2, delta=0.01):
    '''Find a point where seqs1 and seqs2 both have high density and filter
    everything not within delta fraction of that point.

    Over all intervals [low,high], where high = low*(1+delta), finds the
    interval of energies that maximizes the quantity

        min_{i=1,2} |{ seq in seqsi : low <= energy(seq) <= high }|

    Assumes both sequence lists have been sorted by binding energy.'''
    lows = seqs1.wcenergies
    highs = lows*(1+delta)
    idxsLow1 = np.arange(len(lows))
    idxsHigh1 = np.searchsorted(seqs1.wcenergies, highs, side='right')
    idxsLow2 = np.searchsorted(seqs2.wcenergies, lows, side='left')
    idxsHigh2 = np.searchsorted(seqs2.wcenergies, highs, side='right')
    sizes1 = idxsHigh1 - idxsLow1
    sizes2 = idxsHigh2 - idxsLow2
    minSizes = np.minimum(sizes1,sizes2)
    idx = np.argmax(minSizes)
    low = seqs1.wcenergies[idx]
    high = low*(1+delta)
    filteredSeqs1 = seqs1.filter_energy(low, high)
    filteredSeqs2 = seqs2.filter_energy(low, high)
    return (filteredSeqs1, filteredSeqs2)



def pack(arr):
    '''Given 2D arr representing bases A,C,G,T as np.ubyte integers 0,1,2,3,
    (same as DNASeqList's seqarr member variable) return numpy array
    representing them 4 bases per np.ubyte.

    This uses a HUGE amount of memory to vectorize operations to create a
    space-efficient representation.'''
    if PYPY: raise NotImplementedError
    numseqs,seqlen = arr.shape
    numbits = seqlen*2
    bits = np.zeros((numseqs,numbits), dtype=np.ubyte)
    bits[:,::2] = (arr & 0b10) >> 1 # get msb of each base
    bits[:,1::2] = (arr & 0b01)     # get lsb of each base
    packedbits = np.packbits(bits,axis=1)
    return packedbits

def unpack(packedbits,seqlen):
    '''Inverse of pack, where seqlen indicates how many total bases per sequence,
    so if it's not a multiple of 4, the last 2, 4, or 6 bits of the unpacked
    bits should be ignored.'''
    if PYPY: raise NotImplementedError
    bits_unchopped = np.unpackbits(packedbits,axis=1)
    numbits = seqlen*2
    numunpackedbits = bits_unchopped.shape[1]
    numbits_tochop = numunpackedbits % numbits
    bits = bits_unchopped[:,:(-numbits_tochop)]
    msbits = bits[:,0:numbits:2]
    lsbits = bits[:,1:numbits:2]
    arr = (msbits << 1) + lsbits
    return arr
