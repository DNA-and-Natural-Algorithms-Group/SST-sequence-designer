'''
Top-level file for DNA single-stranded tile (SST) sequence designer used in the following publication.
 "Diverse and robust molecular algorithms using reprogrammable DNA self-assembly"
 Woods*, Doty*, Myhrvold, Hui, Zhou, Yin, Winfree. (*Joint first co-authors)
With special thanks to Constantine Evans. 

This program takes as input a set of tile types and outputs SSTs encoding them, given as a parameter
file using the syntax given in, for example  ./input_tilesets/ibc_6bit.py. 

See README file for installation info. 

The glues on each tiles should all have "strength 1" in the aTAM. The program has the following design 
criteria (these, and a number of related constraints, are described in more detail in Suppl Info A of the 
above-mentioned publication, and in the above-mentioned example input file):
1) Glues (aka sticky ends) should have roughly equal binding energy (i.e., be isoenergetic).
2) Tiles should have small internal secondary structure.
3) Interactions between all pairs of tiles should be small, whether or not they have equal (complementary) glues.
4) Minimize interactions of sticky ends that are not complementary, but which might end up "close" due to a tile attachment error.
   (one input domain matching, one input domain mismatching) should have very low binding affinity.
5) Minimize interactions of sticky ends that are not complementary, but which might end up "close" by 
   being co-located on a growing lattice frontier. 
6) No GGGG in a tile and no {C,G}^4 in a sticky end.


Regarding criterion (4), the sticky ends might be close because they are on the
same tile type, but this should be handled by case (2). Therefore, we only
assume y and x are "close" if y* is an "output" sticky ends of some tile type, and
x is an "input" sticky ends, and it is possible for them to get close because the tile
where x is could bind by strength 1 to a tile adjacent to y*'s tile. Here is an
illustration (the # represents a boundary between sticky ends):

          /---------#--------->             tile type t1
          |              b*                /---------#--------->
          |                                |    w
          |              w*                |
          \---------#----------            |    y
/----N----#----E---->                      \---------#----------
|
|
|                                           tile type t2
\----W----#----S-----                      /---------#--------->
          /---------#--------->            |    w
          |              y*                |
          |                                |    x
          |              a*                \---------#----------
          \---------#----------

In this example, tile type t1 matches both available output sticky ends w* and
y*. However, tile type t2 matches w* but not y*. Therefore, because of the
match between t2's w and w*, and the mismatch between x and y*, we require
that x and y* have very low binding affinity.

Regarding criterion (5), we would want y* and w* to have little interaction, and w and x to have little interaction. 
'''

from __future__ import print_function

import dsd
import sst_dsd as sd
import random
import string
import collections
import sys
import itertools
import time
import multiprocessing
from multiprocessing.pool import ThreadPool
import math
import re

global_thread_pool = ThreadPool()

_length10 = 10
_length11 = 11

directions = ('S', 'W', 'N', 'E')
canonical_directions = ('N','W')

ALGO_CONF_GLUES_WEIGHT_EXCESS = 4
ALGO_CONF_GLUES_WEIGHT_EXTEND = 2 # Natural number > 1
TILE_SEC_STRUCT_WEIGHT_EXCESS = 4
TILE_SEC_STRUCT_WEIGHT_EXTEND = 2
TILE_PAIRS_HEURISTIC_WEIGHT_EXCESS = 1
TILE_PAIRS_HEURISTIC_WEIGHT_EXTEND = 1
INPUT_GLUE_PAIRS_EXTEND = 1
OUTPUT_GLUE_PAIRS_EXTEND = 1
LATTICE_BINDING_EXTEND = 1

QUIT_EARLY_OPTIMIZATION = True
THREADED = True
THREADED_TILE_PAIRS = True


def is_canonical(direction):
    return direction in canonical_directions

def opposite(direction):
    if   direction == 'S': return 'N'
    elif direction == 'N': return 'S'
    elif direction == 'E': return 'W'
    elif direction == 'W': return 'E'
    else:                  raise ValueError('%s is not a direction' % direction)

def direction2axis(direction):
    if   direction == 'S': return 'NS'
    elif direction == 'N': return 'NS'
    elif direction == 'E': return 'EW'
    elif direction == 'W': return 'EW'
    else:                  raise ValueError('%s is not a direction' % direction)

# maketrans is in string if python2 but in str if python3
try:
    _wctable = str.maketrans('ACGTacgt','TGCAtgca')
except:
    _wctable = string.maketrans('ACGTacgt','TGCAtgca')

def wc(seq):
    '''Return reverse Watson-Crick complement of seq'''
    return seq.translate(_wctable)[::-1]

def random_iter(iterable):
    "Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    return iter(random.sample(pool, len(pool)))


class Glue:
    _glues = dict()

    def factory(label, axis, end_constraint):
        '''Create new glue if none already exists with given label, axis, and lower
        and upper bounds on strength (lowDG and highDG), length, and biotin_direction direction.'''
        if (label,axis) in Glue._glues:
            glue = Glue._glues[(label,axis)]
            if end_constraint != glue.end_constraint:
                raise ValueError('glue with label %s and axis %s already exists with different EndConstraint' % (label,axis))
        else:
            new_glue = Glue(label, axis, end_constraint)
            Glue._glues[(label,axis)] = new_glue
            glue = new_glue
        return glue

    factory = staticmethod(factory)

    def __init__(self, label, axis, end_constraint):
        '''axis is either "NS" or "EW"; label is a string'''
        self.label = label
        self.end = None
        self.axis = axis
        self.end_constraint = end_constraint
        end_constraint.glues.append(self)
        self.tiles = []
        self.end_bound = False
        self.conflicting_glues_first_order = set()
        self.conflicting_glues_generalized = set()

    def __str__(self):
        return self.label + ';' + self.axis

    def __repr__(self):
        return str(self)

#     def __repr__(self): return 'Glue(label=%s,end=%s,axis=%s)' % (self.label,self.end,self.axis)

    def __hash__(self):
        return hash(self.label + ';AXIS;' + self.axis)

    def __eq__(self, other):
        return self.label == other.label and self.axis == other.axis

    def __cmp__(self, other):
        cmp_label = cmp(self.label, other.label)
        if cmp_label != 0:
            return cmp_label
        else:
            return cmp(self.axis, other.axis)

    def add_tile(self, tile):
        if tile not in self.tiles:
            self.tiles.append(tile)

    def input_tiles(self):
        '''Return tiles for which this is an input glue'''
        return [tile for tile in self.tiles if self in [tile.glues['N'], tile.glues['W']]]

    def output_tiles(self):
        '''Return tiles for which this is an output glue'''
        return [tile for tile in self.tiles if self in [tile.glues['S'], tile.glues['E']]]

    def canonical_tiles(self):
        '''Return tiles for which this is a glue in a canonical direction'''
        return [tile for tile in self.tiles if self in [tile.glues[direction] for direction in canonical_directions]]

    def canonical_direction(self):
        if self.axis == 'NS':
            for direction in ['N','S']:
                if direction in canonical_directions:
                    return direction
        else:
            for direction in ['E','W']:
                if direction in canonical_directions:
                    return direction


    def direction(self, which):
        if   self.axis == 'NS' and which == 0: return 'S'
        elif self.axis == 'NS' and which == 1: return 'N'
        elif self.axis == 'EW' and which == 0: return 'W'
        elif self.axis == 'EW' and which == 1: return 'E'
        else: raise ValueError('which must be 0 or 1, not %s' % which)

    def _check_direction(self, direction):
        if self.axis == 'NS' and direction not in ('N','S'):
            raise ValueError('glue axis %s must be used with direction either N or S, %s is invalid' % (self.axis, direction))
        elif self.axis == 'EW' and direction not in ('E','W'):
            raise ValueError('glue axis %s must be used with direction either E or W, %s is invalid' % (self.axis, direction))

    def get_parity(self, direction):
        '''Get parity of tile that has this glue in given direction.'''
        for tile in self.tiles:
            if tile.glue(direction) == self:
                return tile.parity
        raise ValueError('no tile has glue %s in direction %s' % (repr(self), direction))

    def input_parity(self):
        for tile in self.tiles:
            if tile.glue('N') == self:
                return tile.parity
            elif tile.glue('W') == self:
                return tile.parity
            elif tile.glue('S') == self:
                return 1 - tile.parity
            elif tile.glue('E') == self:
                return 1 - tile.parity
        raise ValueError('no tile has glue %s in any direction' % (repr(self)))

    def get_end(self, direction, include_biotins=False):
        if self.end is None:  return None
        self._check_direction(direction)

        if is_canonical(direction):
            end = self.end
        else:
            end = wc(self.end)

        if include_biotins and self.end_constraint.biotin_direction == direction:
            if not self.tiles:
                raise ValueError('assign this glue to a tile before calling get_end with include_biotins=True')
            parity = self.get_parity(direction)
            biotin_pos = biotin_end_pos(direction, parity)
            if end[biotin_pos] != 'T':
                raise ValueError('to include internal biotin_direction on glue %s of tiles %s, base at position %d must be T, but instead is %s' %
                                 (self, self.tiles[0], biotin_pos, end[biotin_pos]))
            end = end[:biotin_pos] + '/iBiodT/' + end[biotin_pos+1:]
        return end

    def input_end(self):
        '''If axis is "NS", return end of "N", if axis is "EW", return end of "W"'''
        if not self.assigned():
            raise ValueError('glue %s is not assigned' % self.label)
        if self.axis == 'NS':
            return self.get_end('N')
        else:
            return self.get_end('W')

    def output_end(self):
        '''If axis is "NS", return end of "S", if axis is "EW", return end of "E"'''
        if not self.assigned():
            raise ValueError('glue %s is not assigned' % self.label)
        if self.axis == 'NS':
            return self.get_end('S')
        else:
            return self.get_end('E')

    def canonical_end(self):
        '''If axis is "NS", return end of "S", if axis is "EW", return end of "W"'''
        if not self.assigned():
            raise ValueError('glue %s is not assigned' % self.label)
        return self.end

    def bind(self, end):
        '''Bind this glue to given end, so that it cannot be assigned anything else.

        Unlike assign, bind is used to state that a glue MUST have the given end;
        this is used to design new ends for glues that are part of a system
        of pre-existing tiles, where some ends are already existing and cannot
        be changed.'''
        if self.end_bound:
            raise ValueError('glue %s is already bound to end %s; cannot re-bind to end %s') % (self.label, self.end, end)
        self.end_bound = True
        self._assign_no_check_for_bind(end)

    def _assign_no_check_for_bind(self, end):
        if self.tiles:
            tile = self.tiles[0]
            self_direction = None
            for direction in directions:
                if tile.glue(direction) == self:
                    self_direction = direction
                    break
            proper_len = tile.len_end(self_direction)
            if len(end) != proper_len:
                raise ValueError('glue %s cannot be assigned to end %s; length %d is required but end is of length %d' % (self.label, end, proper_len, len(end)))
        self.end = end

    def assign(self, end):
        if self.end_bound:
            raise ValueError('glue %s is bound and cannot be assigned end %s' % (self.label,end))
        self._assign_no_check_for_bind(end)

    def unassign(self):
        self.end = None

    def bound(self):
        return self.end_bound

    def assigned(self):
        return self.end is not None


def short_list(lst):
    if len(lst) < 1000:
        return str(lst)
    else:
        return '%s ... %s' % (str(lst[:10])[:-1], str(lst[-10:])[1:])


class EndConstraint:
    '''This represents a constraint on a glue that affects which sticky ends it
    may have (e.g., glue label does not constrain sticky end, but length does).

    Currently it refers to:
    1) lowDG and highDG, bounds on the energy of the end assigned to the glue
    2) the length end assigned to the glue
    3) (potential) direction of biotin_direction on the glue (which means that base must
       be either a T or an A, depending on whether it is canonical)
    4) the alphabet to use (for three-letter codes)'''
    _end_constraints = dict()

    def factory(lowDG, highDG, length, biotin_direction, endGC, endAT, orth_any, orth_any_ave, orth_algorithmic_conflict, orth_colocated, domain_indv_sec_struct, domain_pair_sec_struct, domain_pair_sec_struct_ave, tile_sec_struct, hamming, orth_algorithmic_conflict_generalized, lattice_binding_lower_threshold): #, seqs):
        '''Create new end constraint if none already exists.'''
        args = tuple((lowDG, highDG, length, biotin_direction, endGC, endAT,
                      orth_any, orth_any_ave, orth_algorithmic_conflict, orth_colocated, domain_indv_sec_struct,
                      domain_pair_sec_struct, domain_pair_sec_struct_ave, tile_sec_struct, hamming, orth_algorithmic_conflict_generalized, lattice_binding_lower_threshold))
        if args in EndConstraint._end_constraints:
            end_constraint = EndConstraint._end_constraints[args]
        else:
            new_end_constraint = EndConstraint(*args)
            EndConstraint._end_constraints[args] = new_end_constraint
            end_constraint = new_end_constraint
        return end_constraint

    factory = staticmethod(factory)

    def end_constraints():
        '''Return list of all EndConstraints created in system.'''
        return list(EndConstraint._end_constraints.values())

    end_constraints = staticmethod(end_constraints)

    def __init__(self, lowDG, highDG, length, biotin_direction, endGC, endAT,
                 orth_any, orth_any_ave, orth_algorithmic_conflict, orth_colocated, domain_indv_sec_struct,
                 domain_pair_sec_struct, domain_pair_sec_struct_ave, tile_sec_struct,
                 hamming, orth_algorithmic_conflict_generalized, lattice_binding_lower_threshold): #, seqs):
        self.lowDG = lowDG
        self.highDG = highDG
        self.lowPF = None
        self.highPF = None
        self.length = length
        self.biotin_direction = biotin_direction
        self.endGC = endGC
        self.endAT = endAT
        self.orth_any = orth_any
        self.orth_any_ave = orth_any_ave
        self.orth_algorithmic_conflict = orth_algorithmic_conflict
        self.orth_colocated = orth_colocated
        self.domain_indv_sec_struct = domain_indv_sec_struct
        self.domain_pair_sec_struct = domain_pair_sec_struct
        self.domain_pair_sec_struct_ave = domain_pair_sec_struct_ave
        self.tile_sec_struct = tile_sec_struct
        self.hamming = hamming
        self.orth_algorithmic_conflict_generalized = orth_algorithmic_conflict_generalized
        self.lattice_binding_lower_threshold = lattice_binding_lower_threshold
#         self.seqs = seqs
        self.glues = []

    def __str__(self):
        return 'EndConstraint(lowDG=%.1f,highDG=%.1f,len=%d,      bio_dir=%s,            endGC=%s,   endAT=%s,   orth_any=%.1f, orth_any_ave=%.1f,orth_algorithmic_conflict=%.1f, orth_colocated=%.1f, domain_indv_sec_struct=%.1f, domain_pair_sec_struct=%.1f, domain_pair_sec_struct_ave=%.1f, tile_sec_struct=%.1f, hamming=%d,   orth_algorithmic_conflict_generalized=%.1f, lattice_binding_lower_threshold=%.1f)' \
            %           (self.lowDG,self.highDG,self.length, self.biotin_direction, self.endGC, self.endAT, self.orth_any, self.orth_any_ave,self.orth_algorithmic_conflict, self.orth_colocated, self.domain_indv_sec_struct, self.domain_pair_sec_struct, self.domain_pair_sec_struct_ave, self.tile_sec_struct, self.hamming, self.orth_algorithmic_conflict_generalized, self.lattice_binding_lower_threshold)

    def __hash__(self):
        if not hasattr(self,'_hash'):
            self._hash = hash(str(self))
        return self._hash

    def __eq__(self, other):
        return self is other
#         return hash(self) == hash(other)

    def __cmp__(self, other):
        return cmp(hash(self),hash(other))

class Tile:
    def __init__(self, name, s, w, n, e, parity, sec_struct, orth, orth_share):
        '''Create tile with given name, glues, and "SST parity".

        name       = name of tile
        s, w, n, e = glues
        parity     = 0 if inner sticky ends are length 10,
                     1 if inner sticky ends are length 11 '''
        self.name = name
        self.glues = {'N': n, 'S': s, 'E': e, 'W': w}
        self.glue_set = set([n,s,e,w])
        for glue in [n,s,e,w]:
            glue.add_tile(self)
        self.parity = parity
        self.sec_struct = sec_struct
        self.orth = orth
        self.orth_share = orth_share
        for (direction,glue) in list(self.glues.items()):
            if glue.bound():
                seq = glue.get_end(direction)
                seq_len = len(seq)
                if not seq_len == len_end(direction, parity):
                    raise ValueError('sequence %s of glue %s is wrong length for direction %s on tile %s with SST parity %d' % (seq, repr(glue), direction, self.name, parity))

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __cmp__(self, other):
        return cmp(self.name, other.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return 'Tile(name=%s,S="%s",W="%s",N="%s",E="%s",parity=%d)' % (self.name,
                                                                self.glues['S'],
                                                                self.glues['W'],
                                                                self.glues['N'],
                                                                self.glues['E'],
                                                                self.parity)

    def idt_sequence(self):
        '''DNA sequence, together with possible biotin_direction labeling.'''
        seq = self.sequence_spaced(include_biotins=True)
        #if self.biotin_direction: seq = '/5Biosg/ ' + seq
        return seq

    def idt(self):
        '''IDT order string for this tile.'''
        return '%s;%s;%s' % (self.name, self.idt_sequence(), 'sst')

    def glue(self, direction):
        if direction not in directions:
            raise ValueError('%s is not a direction' % direction)
        return self.glues[direction]

    def all_glues(self):
        '''Return all glues on this tile.'''
        return list(self.glues.values())

    def unbound_glues(self):
        '''Return all unbound glues on this tile.'''
        uncons_glues = []
        for glue in self.glues.values():
            if not glue.bound():
                uncons_glues.append(glue)
        return uncons_glues

    def len_end(self, direction):
        return len_end(direction, self.parity)

    def directions_with_length(self, length):
        '''Return directions of glues with given length'''
        directions = []
        for direction in directions:
            if self.len_end(direction) == length: directions.append(direction)
        return directions

    def _ends_in_order(self, include_biotins=False):
        '''List of DNA sequence of sticky ends in order S, W, N, E.'''
        ends = []
        for direction in ['S','W','N','E']:
            glue = self.glue(direction)
            end = glue.get_end(direction, include_biotins=include_biotins)
            if not end: end = '?' * self.len_end(direction)
            ends.append(end)
        return ends

    def sequence(self, include_biotins=False):
        '''DNA sequence of sticky ends in order S, W, N, E.'''
        return ''.join(self._ends_in_order(include_biotins))

    def sequence_spaced(self, include_biotins=False):
        '''DNA sequence of sticky ends in order S, W, N, E, with spaces between ends.'''
        return ' '.join(self._ends_in_order(include_biotins))

def len_end(direction, parity):
    if parity == 1:
        if   direction in ['S','E']: return _length10
        elif direction in ['N','W']: return _length11
        else: raise ValueError('%s is not a direction' % direction)
    elif parity == 0:
        if   direction in ['S','E']: return _length11
        elif direction in ['N','W']: return _length10
        else: raise ValueError('%s is not a direction' % direction)

def biotin_tile_pos(direction):
    '''Which position should the internal dT with biotin_direction go, if the end lies in
    the given direction. (Should make the major groove, where the internal
    biotin_direction resides,  point towards the "inside" of the SST tube).

    The indexing is with respect to the entire tile, not with respect to the
    start of the glue end (hence the positions range between 0 and 41, not
    between 0 and 9 or 0 and 10).'''
    if _length10 == 10 and _length11 == 11:
        if   direction == 'S': return 5
        elif direction == 'W': return 15
        elif direction == 'N': return 27
        elif direction == 'E': return 37
        else: raise ValueError('%s is not a valid direction' % direction)
    else: # this is for debugging purposes when I want to test fewer sticky ends than 4^10 or 4^11
        if   direction == 'S': return 2
        elif direction == 'W': return 8
        elif direction == 'N': return 13
        elif direction == 'E': return 18
        else: raise ValueError('%s is not a valid direction' % direction)

def biotin_end_pos(direction, parity):
    '''Like biotin_tile_pos, but relative to the glue end, not the
    whole tile, so indexed from 0 to 9 or 10. Depends on "SST parity" of tile
    since that determines the length of the glue end.

    Note that if the end is not canonical, to have a T at a
    position pos, needs its canonical partner to have an A at position
    end_length - pos.

    parity     = 0 if inner sticky ends ('W','N') are length 10,
                 1 if inner sticky ends are length 11'''
    tile_pos = biotin_tile_pos(direction)
    if   direction == 'S': return tile_pos
    elif direction == 'W': return tile_pos - (len_end('S', parity))
    elif direction == 'N': return tile_pos - (len_end('S', parity) + len_end('W', parity))
    elif direction == 'E': return tile_pos - (len_end('S', parity) + len_end('W', parity) + len_end('N', parity))

def random_bad_unconstrained_tile(bad_tiles):
    bad_unconstrained_tiles = []
    for tile in bad_tiles:
        if tile.unbound_glues():
            bad_unconstrained_tiles.append(tile)
    return random.choice(bad_unconstrained_tiles)


def split(l, num):
    """Yield successive chunks from l, each of same size, to make num total.

    http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python"""
    assert 0 < num <= len(l)
    n = len(l) // num
    if len(l) % num != 0:
        n += 1
    for i in range(0, len(l), n):
        yield l[i:i+n]

LOG_BAD_GLUES = False
def log_bad_glues(reason):
    if LOG_BAD_GLUES:
        print(reason + ';', end=' ')

class TileSet:
    def __init__(self, tiles, detect_nondeterminism, mutually_exclusive_tilename_pairs):
        for tile in tiles:
            tile.tileset = self
        self.tiles = tiles
        tile_name_list = [tile.name for tile in tiles]
        tile_name_set = set(tile_name_list)
        if len(tile_name_set) != len(tile_name_list):
            raise ValueError('duplicate tile name')
        self.name2tile = { tile.name: tile for tile in tiles }
        self.tile_pairs_to_check = [(t1,t2) for (t1,t2) in itertools.combinations_with_replacement(tiles,2)
                                    if (t1.name,t2.name) not in mutually_exclusive_tilename_pairs]
        self.tile_pair_names_to_check = [ (t1.name,t2.name) for (t1,t2) in self.tile_pairs_to_check ]
        self.find_conflicting_glues(detect_nondeterminism)
        self.check_parity_against_glues()

    def end_constraints_unbound(self):
        if not hasattr(self,'_end_constraints_unbound'):
            set_end_constraints_unbound = set()
            for glue in self.glues_unbound():
                set_end_constraints_unbound.add(glue.end_constraint)
            self._end_constraints_unbound = list(set_end_constraints_unbound)
        return self._end_constraints_unbound

    def check_parity_against_glues(self):
        '''Check that the tiles are consistent in that glues only join tiles
        to tiles of the opposite parity.'''
        for tile in self.tiles:
            for direc in directions:
                glue = tile.glue(direc)
                opp = opposite(direc)
                for other_tile in glue.tiles:
                    if other_tile.glue(opp) == glue and tile.parity == other_tile.parity:
                        raise ValueError('glue %s joins tiles %s and %s but they have the same parity' % (glue,tile,other_tile))

    def find_conflicting_glues(self, detect_nondeterminism):
        '''Search through all tiles, assuming north and west are "input" sticky
        ends, and assign to each such glue a list of its conflicting glues.
        These are glues it could end up next to because
        1) the other input glue is shared with another tile,
        2) two output glues of a given tile,
        3) two output glues of diagonally adjacent tiles (e.g., w* and y* above).'''
        # reason 1
        for (tile1,tile2) in self.tile_pairs_to_check:
            if tile1 == tile2:
                continue
            ng1 = tile1.glue('N')
            wg1 = tile1.glue('W')
            ng2 = tile2.glue('N')
            wg2 = tile2.glue('W')
            if detect_nondeterminism and ng1 == ng2 and wg1 == wg2:
                raise ValueError('nondeterminism detected: tiles %s and %s share north glue %s and west glue %s in common' % (tile1.name,tile2.name,ng1.label,ng2.label))

            # In the following case we have a nondeterministic tileset, that is:
            # tile1 != tile2, and tile1 & tile2
            # share both N and W glues (i.e. tile1 & tile2 share both input sides)
            if ng1 == ng2 and wg1 == wg2:
                continue
            if ng1 == ng2:
                wg1.conflicting_glues_first_order.add(wg2)
                wg2.conflicting_glues_first_order.add(wg1)
            if wg1 == wg2:
                ng1.conflicting_glues_first_order.add(ng2)
                ng2.conflicting_glues_first_order.add(ng1)


            if len(tile1.name.split(';')) == 3 and len(tile2.name.split(';')) == 3:
                # If tile names are of the form ___;___;___ then the code assumes they are using
                # our proofreading block naming convention and carries out an orthogonality test
                # for generalised algorithmic (tile attachment) errors. Generalised algorithmic (tile attachment) 
                # errors are errors where we assume one or more errors have already occurred elsewhere
                # in a proof-reading block. Generalised algorithmic (tile attachment) errors are also described in Suppl. Info. A.  
                gate_position1,_,proof_internal_pos1 = tile1.name.split(';')
                gate_position2,_,proof_internal_pos2 = tile2.name.split(';')
                if gate_position1 == gate_position2 and proof_internal_pos1 == proof_internal_pos2:
                    if wg1 != wg2:
                        wg1.conflicting_glues_generalized.add(wg2)
                        wg2.conflicting_glues_generalized.add(wg1)
                    if ng1 != ng2:
                        ng1.conflicting_glues_generalized.add(ng2)
                        ng2.conflicting_glues_generalized.add(ng1)



    def check4G(self):
        # find glues whose tiles have 4 G's
        bad_glues = []
        for tile in self.tiles:
            if tile.unbound_glues(): # if all bound, then nothing we can do
                seq = tile.sequence()
                if 'GGGG' in seq:
                    log_bad_glues('4G')
                    for dir1, dir2 in [('S', 'W'), ('W', 'N'), ('N', 'E')]:
                        g1 = tile.glue(dir1)
                        g2 = tile.glue(dir2)
                        end1 = g1.get_end(dir1)
                        end2 = g2.get_end(dir2)
                        if 'GGGG' in end1 + end2:
                            if not g1.bound():
                                bad_glues.append(g1)
                            if not g2.bound():
                                bad_glues.append(g2)
                            if g1.bound() and g2.bound():
                                print('\n' + '*' * 79 + 'WARNING: glues %s and %s on tile %s are bound to sequences %s and %s, which produces a G quadruplex\n' % (g1, g2, tile, end1, end2))
        return bad_glues

    def check_every_glue(self, temperature):
        bad_glues = []
        # check all unbound glues for too much interaction with themselves
        #   (both its end and the WC complement)
        for glue_unbound in self.glues_unbound():
            ec = glue_unbound.end_constraint
            end_unbound = glue_unbound.input_end()
            end_unbound_wc = glue_unbound.output_end()
            energy = sd.binding(end_unbound, end_unbound, temperature)
            if energy > ec.orth_any:
                log_bad_glues('oa%.1f' % energy)
                bad_glues.append(glue_unbound)
            energy = sd.binding(end_unbound_wc, end_unbound_wc, temperature)
            if energy > ec.orth_any:
                log_bad_glues('oa%.1f' % energy)
                bad_glues.append(glue_unbound) # check this unbound glue against all glues (potentially bound)
            # for too much interaction; if found, add unbound glue (also add glue
            #   if it is not bound)
            for glue in self.glues():
                end = glue.input_end()
                end_wc = glue.output_end()
                problem = False
                if end != end_unbound_wc:
                    energy = sd.binding(end, end_unbound, temperature)
                    if energy > ec.orth_any:
                        problem = True
                        log_bad_glues('oa%.1f' % energy)
                    energy = sd.binding(end_wc, end_unbound_wc, temperature)
                    if energy > ec.orth_any:
                        problem = True
                        log_bad_glues('oa%.1f' % energy)
                if end != end_unbound:
                    energy = sd.binding(end, end_unbound_wc, temperature)
                    if energy > ec.orth_any:
                        problem = True
                        log_bad_glues('oa%.1f' % energy)
                    energy = sd.binding(end_wc, end_unbound, temperature)
                    if energy > ec.orth_any:
                        problem = True
                        log_bad_glues('oa%.1f' % energy)
                if problem:
                    bad_glues.append(glue_unbound)
                    if not glue.bound():
                        bad_glues.append(glue)

        return bad_glues

    def check_tile_sec_struct(self, temperature, num_bad_now, num_bad_opt, threaded, weight_excess):
        '''Find glues whose tiles have too much secondary structure.'''
        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
            bad_glues = []
            if threaded:
                group_size = global_thread_pool._processes
                for tile_group in grouper(self.tiles, group_size):
                    tile_group = [tile for tile in tile_group if tile] # strip out None's that fill in to make size of last group equal to group_size
                    results = [global_thread_pool.apply_async(sd.hairpin, args=(tile.sequence(), temperature)) for tile in tile_group]
                    energies = [result.get() for result in results]
                    for tile,energy in zip(tile_group, energies):
                        if tile.unbound_glues() and tile.sec_struct > 0: # if all bound, then nothing we can do
                            if energy > tile.sec_struct:
                                # add proportional to excess of energy
                                num_to_add = int(math.ceil(weight_excess*(energy - tile.sec_struct)))
                                bad_glues.extend(tile.unbound_glues() * num_to_add)
                                if len(bad_glues)*TILE_SEC_STRUCT_WEIGHT_EXTEND + num_bad_now > num_bad_opt:
                                    return bad_glues
            else:
                for tile in self.tiles:
                    energy = sd.hairpin(tile.sequence(), temperature) if tile.sec_struct > 0 else 0
                    if tile.unbound_glues() and tile.sec_struct > 0: # if all bound, then nothing we can do
                        if energy > tile.sec_struct:
                            # add proportional to excess of energy
                            num_to_add = int(math.ceil(weight_excess*(energy - tile.sec_struct)))
                            bad_glues.extend(tile.unbound_glues() * num_to_add)
                            if len(bad_glues)*TILE_SEC_STRUCT_WEIGHT_EXTEND + num_bad_now > num_bad_opt:
                                return bad_glues
            return bad_glues
        else:
            if threaded:
                def ss_if_positive(tile, temperature):
                    if tile.sec_struct < 0:
                        return (tile,0)
                    else:
                        return (tile,sd.hairpin(tile.sequence(), temperature))
                results = [global_thread_pool.apply_async(ss_if_positive, args=(tile, temperature)) for tile in self.tiles]
                tile_energies = [result.get() for result in results]
            else:
                tile_energies = []
#                 tile_num = 1
                for tile in self.tiles:
#                     print 'checking tile# {}'.format(tile_num)
#                     tile_num += 1
                    energy = sd.hairpin(tile.sequence(), temperature) if tile.sec_struct > 0 else 0
                    tile_energies.append((tile,energy))

            bad_glues = []
            for tile,energy in tile_energies:
                if tile.unbound_glues() and tile.sec_struct > 0: # if all bound, then nothing we can do
                    if energy > tile.sec_struct:
                        # add proportional to excess of energy
                        num_to_add = int(math.ceil(weight_excess*(energy - tile.sec_struct)))
                        bad_glues.extend(tile.unbound_glues() * num_to_add)
            return bad_glues

    def check_input_pairs(self, temperature, num_bad_now, num_bad_opt, threaded):
        '''Check pairs of ends that each tile's input ends bind to (which must
        be pulled apart for the tile to bind).'''
        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
            bad_glues = []
            for tile in self.tiles:
                ng = tile.glue('N')
                wg = tile.glue('W')
                if ng.bound() and wg.bound():
                    continue
                end1 = ng.output_end()
                end2 = wg.output_end()
                orth_colocated = ng.end_constraint.orth_colocated
                if orth_colocated < 0:
                    continue
                energy = eval_colocated_end_pair(end1, end2, temperature)
                if energy > orth_colocated:
                    if not ng.bound(): bad_glues.append(ng)
                    if not wg.bound(): bad_glues.append(wg)
                    if len(bad_glues)*INPUT_GLUE_PAIRS_EXTEND + num_bad_now > num_bad_opt:
                        return bad_glues
            return bad_glues
        else:
            glue_end_pairs = []
            for tile in self.tiles:
                ng = tile.glue('N')
                wg = tile.glue('W')
                if ng.bound() and wg.bound():
                    continue
                nend = ng.output_end()
                wend = wg.output_end()
                orth_colocated = ng.end_constraint.orth_colocated
                if orth_colocated < 0:
                    continue
                glue_end_pairs.append((ng,wg,nend,wend,orth_colocated))
            end_pairs = [(nend,wend) for (_,_,nend,wend,orth_colocated) in glue_end_pairs]
            energies = eval_colocated_end_pairs(end_pairs, temperature, threaded)
            bad_glues = []
            for ((ng, wg, nend, wend, orth_colocated), energy) in zip(glue_end_pairs, energies):
                if energy > orth_colocated:
                    if not ng.bound(): bad_glues.append(ng)
                    if not wg.bound(): bad_glues.append(wg)
            return bad_glues

    def check_output_pairs(self, temperature, num_bad_now, num_bad_opt, threaded):
        '''Check pairs of ends that are output ends of tiles (which must
        be pulled apart another tile to bind using either of them).'''
        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
            bad_glues = []
            for tile in self.tiles:
                sg = tile.glue('S')
                eg = tile.glue('E')
                if sg.bound() and eg.bound():
                    continue
                end1 = sg.output_end()
                end2 = eg.output_end()
                orth_colocated = sg.end_constraint.orth_colocated
                if orth_colocated < 0:
                    continue
                energy = eval_colocated_end_pair(end1, end2, temperature)
                if energy > orth_colocated:
                    if not sg.bound(): bad_glues.append(sg)
                    if not eg.bound(): bad_glues.append(eg)
                    if len(bad_glues)*OUTPUT_GLUE_PAIRS_EXTEND + num_bad_now > num_bad_opt:
                        return bad_glues
            return bad_glues
        else:
            glue_end_pairs = []
            for tile in self.tiles:
                sg = tile.glue('S')
                eg = tile.glue('E')
                if sg.bound() and eg.bound():
                    continue
                out_end1 = sg.output_end()
                out_end2 = eg.output_end()
                orth_colocated = sg.end_constraint.orth_colocated
                if orth_colocated < 0:
                    continue
                glue_end_pairs.append((sg,eg,out_end1,out_end2,orth_colocated))

    #             energy = sd.binding(out_end1, out_end2, temperature)
            end_pairs = [(end1,end2) for (_,_,end1,end2,orth_colocated) in glue_end_pairs]
            energies = eval_colocated_end_pairs(end_pairs, temperature, threaded)
            bad_glues = []
            for ((sg, eg, end1, end2, orth_colocated), energy) in zip(glue_end_pairs, energies):
                if energy > orth_colocated:
                    if not sg.bound(): bad_glues.append(sg)
                    if not eg.bound(): bad_glues.append(eg)
            return bad_glues

    def check_lattice_binding(self, temperature, num_bad_now, num_bad_opt, threaded):
        '''Check the strnegth of binding of each tile, by its two input ends,
        to the matching pair lattice ends.'''
        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
            bad_glues = []
            for tile in self.tiles:
                wg = tile.glue('W')
                ng = tile.glue('N')
                if wg.bound() and ng.bound():
                    continue
                end1 = wg.input_end() # was "output_end"
                end2 = ng.input_end()
                lower_threshold_binding_energy = ng.end_constraint.lattice_binding_lower_threshold # assume all glues get same lattice_binding_lower_threshold
                if lower_threshold_binding_energy < 0:
                    continue

                energy = eval_lattice_binding_energy(end1, end2, temperature)
                if energy < lower_threshold_binding_energy:
                    if not wg.bound(): bad_glues.append(wg)
                    if not ng.bound(): bad_glues.append(ng)
                    if len(bad_glues)*LATTICE_BINDING_EXTEND + num_bad_now > num_bad_opt:
                        return bad_glues
            return bad_glues
        else:
            glue_end_pairs = []
            for tile in self.tiles:
                wg = tile.glue('W')
                ng = tile.glue('N')
                if wg.bound() and ng.bound():
                    continue
                end1 = wg.input_end()
                end2 = ng.input_end()
                lower_threshold_binding_energy = ng.end_constraint.lattice_binding_lower_threshold # assume all glues get same lattice_binding_lower_threshold
                if lower_threshold_binding_energy < 0:
                    continue
                glue_end_pairs.append((wg,ng,end1,end2,lower_threshold_binding_energy))

            end_pairs = [(end1,end2) for (_,_,end1,end2,lower_threshold_binding_energy) in glue_end_pairs]
            energies = eval_lattice_binding_energies(end_pairs, temperature, threaded)
            bad_glues = []
            for ((wg, ng, end1, end2, lower_threshold_binding_energy), energy) in zip(glue_end_pairs, energies):
                if energy < lower_threshold_binding_energy:
                    if not wg.bound(): bad_glues.append(wg)
                    if not ng.bound(): bad_glues.append(ng)
            return bad_glues

    def check_algorithmic_conflicting_glues(self, temperature, num_bad_now, num_bad_opt, threaded, weight_excess):
        ''' Find glues that have too much interaction with an algorithmic conflicting glue

        (could end up near each other during a strength-1 binding event).'''
        bad_glues_first_order = []
        bad_glues_generalized = []
        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
#             bad_glues = []
            evals = 1
            for glue in self.glues():
                end = glue.input_end()
                for conflicting_glue in glue.conflicting_glues_first_order:
                    conflicting_end = conflicting_glue.output_end()
                    orth_algorithmic_conflict_first_order = conflicting_glue.end_constraint.orth_algorithmic_conflict
                    if orth_algorithmic_conflict_first_order < 0:
                        continue
                    # the proper order in which to concatenate the ends depends
                    # on the axis, observe:
                    #           /---------#--------->
                    #           |
                    #           |
                    #           |
                    #           \----W----#----------
                    # /---------#----E----> this should be W (in) then E (out)
                    # |
                    # |
                    # |
                    # \---------#----S----- this should be S (out) then N (in)
                    #           /----N----#--------->
                    #           |
                    #           |
                    #           |
                    #           \---------#----------

                    if glue.axis == 'EW': end1,end2 = end, conflicting_end
                    else:                 end1,end2 = conflicting_end, end

                    energy = eval_colocated_end_pair(end1, end2, temperature)
                    evals += 1
                    if energy > orth_algorithmic_conflict_first_order:
                        num_to_add = int(math.ceil(weight_excess*(energy - orth_algorithmic_conflict_first_order)))
                        if not glue.bound():
                            bad_glues_first_order.extend([glue] * num_to_add)
                        if not conflicting_glue.bound():
                            bad_glues_first_order.extend([conflicting_glue] * num_to_add)
#                         print 'testing for early quit; len(bad_glues) = {}; num_bad_now = {}; len(bad_glues)+num_bad_now = {}; num_bad_opt = {}'.format(len(bad_glues), num_bad_now, len(bad_glues)+num_bad_now, num_bad_opt)
                        if len(bad_glues_first_order)*ALGO_CONF_GLUES_WEIGHT_EXTEND + num_bad_now > num_bad_opt:
#                             print '\nnum acg evaluations: {}'.format(evals)
                            return (bad_glues_first_order, bad_glues_generalized) # bad_glues_first_order
                for conflicting_glue in glue.conflicting_glues_generalized:
                    conflicting_end = conflicting_glue.output_end()
                    orth_algorithmic_conflict_generalized = conflicting_glue.end_constraint.orth_algorithmic_conflict_generalized
                    if orth_algorithmic_conflict_generalized < 0:
                        continue

                    if glue.axis == 'EW': end1,end2 = end, conflicting_end
                    else:                 end1,end2 = conflicting_end, end

                    energy = eval_colocated_end_pair(end1, end2, temperature)
                    evals += 1
                    if energy > orth_algorithmic_conflict_generalized:
                        num_to_add = int(math.ceil(weight_excess*(energy - orth_algorithmic_conflict_generalized)))
                        if not glue.bound():
                            bad_glues_generalized.extend([glue] * num_to_add)
                        if not conflicting_glue.bound():
                            bad_glues_generalized.extend([conflicting_glue] * num_to_add)
#                         print 'testing for early quit; len(bad_glues) = {}; num_bad_now = {}; len(bad_glues)+num_bad_now = {}; num_bad_opt = {}'.format(len(bad_glues), num_bad_now, len(bad_glues)+num_bad_now, num_bad_opt)
                        if len(bad_glues_generalized)*ALGO_CONF_GLUES_WEIGHT_EXTEND + num_bad_now > num_bad_opt:
#                             print '\nnum acg evaluations: {}'.format(evals)
                            return (bad_glues_first_order, bad_glues_generalized) # bad_glues_generalized
#             print '\nnum acg evaluations: {}'.format(evals)

        else:
            glue_end_pairs_first_order = []
            glue_end_pairs_generalized = []
            for glue in self.glues():
                end = glue.input_end()
                for conflicting_glue in glue.conflicting_glues_first_order:
                    conflicting_end = conflicting_glue.output_end()
                    orth_algorithmic_conflict_first_order = conflicting_glue.end_constraint.orth_algorithmic_conflict
                    if orth_algorithmic_conflict_first_order < 0:
                        continue
                    if glue.axis == 'EW': end1,end2,glue1,glue2 = end, conflicting_end, glue, conflicting_glue
                    else:                 end1,end2,glue1,glue2 = conflicting_end, end, conflicting_glue, glue
                    glue_end_pairs_first_order.append((glue1, glue2, end1, end2, orth_algorithmic_conflict_first_order))
                for conflicting_glue in glue.conflicting_glues_generalized:
                    conflicting_end = conflicting_glue.output_end()
                    orth_algorithmic_conflict_generalized = conflicting_glue.end_constraint.orth_algorithmic_conflict_generalized
                    if orth_algorithmic_conflict_generalized < 0:
                        continue
                    if glue.axis == 'EW': end1,end2,glue1,glue2 = end, conflicting_end, glue, conflicting_glue
                    else:                 end1,end2,glue1,glue2 = conflicting_end, end, conflicting_glue, glue
                    glue_end_pairs_generalized.append((glue1, glue2, end1, end2, orth_algorithmic_conflict_generalized))
            end_pairs_first_order = [(end1,end2) for (_,_,end1,end2,orth_algorithmic_conflict) in glue_end_pairs_first_order]
            end_pairs_generalized = [(end1,end2) for (_,_,end1,end2,orth_algorithmic_conflict) in glue_end_pairs_generalized]
            energies_first_order = eval_colocated_end_pairs(end_pairs_first_order, temperature, threaded)
            energies_generalized = eval_colocated_end_pairs(end_pairs_generalized, temperature, threaded)
            for ((glue1, glue2, end1, end2, orth_algorithmic_conflict), energy) in zip(glue_end_pairs_first_order, energies_first_order):
                if energy > orth_algorithmic_conflict:
                    num_to_add = int(math.ceil(weight_excess*(energy - orth_algorithmic_conflict)))
                    if not glue1.bound():
                        bad_glues_first_order.extend([glue1] * num_to_add)
                    if not glue2.bound():
                        bad_glues_first_order.extend([glue2] * num_to_add)
            for ((glue1, glue2, end1, end2, orth_algorithmic_conflict), energy) in zip(glue_end_pairs_generalized, energies_generalized):
                if energy > orth_algorithmic_conflict:
                    num_to_add = int(math.ceil(weight_excess*(energy - orth_algorithmic_conflict)))
                    if not glue1.bound():
                        bad_glues_generalized.extend([glue1] * num_to_add)
                    if not glue2.bound():
                        bad_glues_generalized.extend([glue2] * num_to_add)

        return (bad_glues_first_order, bad_glues_generalized)

    def populate_tile_pair_name2energy_most_recent(self, temperature, tile_pair_names_to_check_this_time, start_pos, threaded_tile_pairs):
        # DD: this number is chosen because in timing tests, the overhead of
        # RNAduplex is about 18ms, and after that it takes about 1ms per pair
        # of length 42 seqs. So 100 pairs is where overhead of calling
        # RNAduplex is about 20%, i.e., not too large, but a small enough
        # number of pairs to make a bigger payoff for quitting early
        num_pairs_per_call = 100
        num_parallel_lists = global_thread_pool._processes if threaded_tile_pairs else 1
        namepairs = [(n1,n2) for (n1,n2) in tile_pair_names_to_check_this_time[start_pos:start_pos+num_pairs_per_call*num_parallel_lists]]
        seqpairs = [ (self.name2tile[n1].sequence(), self.name2tile[n2].sequence()) for (n1,n2) in namepairs ]

        if threaded_tile_pairs and len(seqpairs) >= num_parallel_lists:
            list_of_list_of_seqpair = list(split(seqpairs, num_parallel_lists))
            lengths = [len(list_of_seqpair) for list_of_seqpair in list_of_list_of_seqpair]
            if sum(lengths) != len(seqpairs):
                raise ValueError('wrong number of tile pairs checked, original list size: {}, sum of lengths of parallel lists: {}'.format(len(seqpairs), sum(lengths)))
            results = [global_thread_pool.apply_async(sd.RNAduplex_multiple, args=(seqpairs, temperature)) for seqpairs in list_of_list_of_seqpair]
            list_of_list_of_energies = [result.get() for result in results]
            energies = itertools.chain.from_iterable(list_of_list_of_energies)
        else:
            energies = sd.RNAduplex_multiple(seqpairs, temperature)

        for ((n1,n2),energy) in zip(namepairs,energies):
            self.tile_pair_name2energy_most_recent[(n1,n2)] = energy
        return start_pos+len(namepairs)

    def check_tile_pairs_heuristic(self, temperature, num_bad_now, num_bad_opt, threaded_tile_pairs, weight_excess):
        if not hasattr(self, 'tile_names_changed_from_opt'): # DW: True iff this is the first call to the function check_tile_pairs_heuristic
            # DW: following will be a dict with keys which are all pairs (i.e. of size 63,190)
            tile_pair_names_to_check_this_time = self.tile_pair_names_to_check
        else:
            tile_pair_names_to_check_this_time = [ (t1name,t2name)
                for (t1name,t2name) in self.tile_pair_names_to_check   # DW: tile_pair_names_to_check is everyone!
                if (t1name in self.tile_names_changed_from_opt
                or  t2name in self.tile_names_changed_from_opt) ]   # DW: tile_names_changed_from_opt was set
                                                                    # in local_search_end_assign(...)
                                                                    # and is the set of tiles (names) that have
                                                                    # changed from the previous best (opt) glue seq set


        if QUIT_EARLY_OPTIMIZATION and num_bad_opt > 0:
            set_tile_pair_names_to_check_this_time = set(tile_pair_names_to_check_this_time)
            self.tile_pair_name2energy_most_recent = dict()
            start_pos = 0
            bad_glues = []
            # assign bad glues based on bad tile pairs
            for (t1name,t2name) in self.tile_pair_names_to_check:
                # if we re-calculated the energy of this tile pair this call, use that, otherwise use the cached value
                if (t1name,t2name) in set_tile_pair_names_to_check_this_time and (t1name,t2name) not in self.tile_pair_name2energy_most_recent:
                    start_pos = self.populate_tile_pair_name2energy_most_recent(temperature, tile_pair_names_to_check_this_time, start_pos, threaded_tile_pairs)
                energy = self.tile_pair_name2energy_most_recent.get((t1name,t2name))
                if energy is None:
                    energy = self.tile_pair_name2energy_opt[(t1name,t2name)]

                t1 = self.name2tile[t1name]
                t2 = self.name2tile[t2name]
                t1orth,t2orth = (t1.orth_share,t2.orth_share) if share_complementary_glues(t1, t2) else (t1.orth,t2.orth)
                if t1orth <= 0 and t2orth <= 0:
                    continue
                if energy > t1orth or energy > t2orth:
                    # add proportional to excess of energy
                    num_to_add = int(math.ceil(weight_excess*(energy - min(t1.orth, t2.orth))))
                    bad_glues.extend((t1.unbound_glues() + t2.unbound_glues()) * num_to_add)
                    if len(bad_glues)*TILE_PAIRS_HEURISTIC_WEIGHT_EXTEND + num_bad_now > num_bad_opt:
                        return bad_glues
            return bad_glues
        else:
            tile_pair_seqs = [(self.name2tile[t1name].sequence(), self.name2tile[t2name].sequence())
                              for (t1name,t2name) in tile_pair_names_to_check_this_time]

            num_parallel_lists = global_thread_pool._processes
            if threaded_tile_pairs and len(tile_pair_seqs) >= num_parallel_lists:
    #             num_parallel_lists = max(global_thread_pool._processes, len(tile_pair_seqs) // 500)
                list_of_list_of_seqpair = list(split(tile_pair_seqs, num_parallel_lists))
                lengths = [len(list_of_seqpair) for list_of_seqpair in list_of_list_of_seqpair]
                if sum(lengths) != len(tile_pair_seqs):
                    raise ValueError('wrong number of tile pairs checked, original list size: {}, sum of lengths of parallel lists: {}'.format(len(tile_pair_seqs), sum(lengths)))
    #             print 'total tile pairs to check: {}\nnum parallel lists: {}\nmax length of each list for each parallel call to {}: {}'.format(len(tile_pair_seqs), num_parallel_lists, 'RNAduplex' if not USE_NUPACK else 'pfunc', lengths[0])
                results = [global_thread_pool.apply_async(sd.RNAduplex_multiple, args=(seqpairs, temperature)) for seqpairs in list_of_list_of_seqpair]
                list_of_list_of_energies = [result.get() for result in results]
                energies = itertools.chain.from_iterable(list_of_list_of_energies)
            else:
                energies = sd.RNAduplex_multiple(tile_pair_seqs, temperature)

            self.tile_pair_name2energy_most_recent = dict(list(zip(tile_pair_names_to_check_this_time, energies)))

            bad_glues = []
            # assign bad glues based on bad tile pairs
    #         start = time.time()
            for (t1name,t2name) in self.tile_pair_names_to_check:
                # if we re-calculated the energy of this tile pair this call, use that, otherwise use the cached value
                energy = self.tile_pair_name2energy_most_recent.get((t1name,t2name))
                if energy is None:
                    energy = self.tile_pair_name2energy_opt[(t1name,t2name)]

                t1 = self.name2tile[t1name]
                t2 = self.name2tile[t2name]
                t1orth,t2orth = (t1.orth_share,t2.orth_share) if share_complementary_glues(t1, t2) else (t1.orth,t2.orth)
                if t1orth <= 0 and t2orth <= 0:
                    continue
                if energy > t1orth or energy > t2orth:
                    # add proportional to excess of energy
                    num_to_add = int(math.ceil(weight_excess*(energy - min(t1.orth, t2.orth))))
                    bad_glues.extend((t1.unbound_glues() + t2.unbound_glues()) * num_to_add)
    #         end = time.time()
    #         print '\n' + '*'*79 + ('\n** %.3f seconds for loop to log bad tile pair glues\n' % ((end - start))) + '*'*79 + '\n'

            return bad_glues

    def check_tile_pairs_heuristic_check_all(self, temperature, threaded_tile_pairs, tile_pair_names_to_check_this_time):
        '''Non-caching version of check_tile_pairs_heuristic for debugging'''
        tile_pair_seqs = [(self.name2tile[t1name].sequence(), self.name2tile[t2name].sequence())
                          for (t1name,t2name) in tile_pair_names_to_check_this_time]

        if threaded_tile_pairs:
            list_of_list_of_seqpair = list(split(tile_pair_seqs, global_thread_pool._processes))
            results = [global_thread_pool.apply_async(sd.RNAduplex_multiple, args=(seqpairs, temperature)) for seqpairs in list_of_list_of_seqpair]
            list_of_list_of_energies = [result.get() for result in results]
            energies = itertools.chain.from_iterable(list_of_list_of_energies)
        else:
            energies = sd.RNAduplex_multiple(tile_pair_seqs, temperature)

        tile_pair_name2energy = dict(list(zip(tile_pair_names_to_check_this_time, energies)))

        return tile_pair_name2energy

        # this code run else; helps me remember variable names
        #self.bad_tile_pairs_opt = self.bad_tile_pairs
        #self.tiles_changed_from_opt = glue.tiles
        #self.name2tile = { tile.name: tile for tile in tiles }
        #self.tile_pairs_to_check = [(t1,t2) for (t1,t2) in itertools.combinations_with_replacement(tiles,2)
        #                            if (t1.name,t2.name) not in mutually_exclusive_tilename_pairs]
        #self.tile_pair_names_to_check = set( (t1.name,t2.name) for (t1,t2) in self.tile_pairs_to_check )


    def log_bad_glue_source(self, num_bad_glues, num_bad_glues_opt, num_bad_glues_orig, last_char):
        width = self.width
        small_width = (width-3)//3
        sys.stdout.write('{:>{width}}{}'.format('{:{small_width}d}#{:{small_width}d}/{:{small_width}d}'.format(num_bad_glues_opt, num_bad_glues, num_bad_glues_orig, small_width=small_width), last_char, width=width))
        sys.stdout.flush()
#         print '\n'*LOG_BAD_GLUES + '# new bad glues due to {:28s} {:5d} / {:5d} originally'.format(src_str, num_bad_glues, num_bad_glues_orig)

    def find_bad_glues(self, temperature, num_bad_opt, threaded=True, threaded_tile_pairs=True, check_tile_pairs=True):
        '''Check all tiles and glues to look for glues that are causing
        unwanted interactions.'''
        global global_thread_pool
        if multiprocessing.cpu_count() != global_thread_pool._processes:
            print('number of cpus changed from %d to %d; allocating new thread pool' % (global_thread_pool._processes, multiprocessing.cpu_count()))
            global_thread_pool = ThreadPool()
        new_bad_glues = []
        self.width = 18
        width=self.width
        print('{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}|{:^{width}}'.format('tile sec struct', 'alg conf 1st', 'alg conf gen', 'lattice binding','in glue pair', 'out glue pair', 'tile pair', 'total', width=width))
        sys.stdout.flush()
#         new_bad_glues.extend(self.check4G())
#         if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt: return new_bad_glues

        CHECK_ONLY_ACG = False

        if not CHECK_ONLY_ACG:
            bad_glues_tss = self.check_tile_sec_struct(temperature, len(new_bad_glues), num_bad_opt, threaded, weight_excess=TILE_SEC_STRUCT_WEIGHT_EXCESS)
            bas_glues_tss = bad_glues_tss*TILE_SEC_STRUCT_WEIGHT_EXTEND
            if not hasattr(self, 'num_bg_tss_orig'): self.num_bg_tss_orig = len(bad_glues_tss)
            if not hasattr(self, 'num_bg_tss_opt'): self.num_bg_tss_opt = -1
            self.num_bg_tss_most_recent = len(bad_glues_tss)
            self.log_bad_glue_source(len(bad_glues_tss), self.num_bg_tss_opt, self.num_bg_tss_orig, '|')
            new_bad_glues.extend(bas_glues_tss)

            if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt:
                self.log_bad_glue_source(-1, self.num_bg_acg_opt_fo, self.num_bg_acg_orig_fo, '|')
                self.log_bad_glue_source(-1, self.num_bg_acg_opt_g, self.num_bg_acg_orig_g, '|')
                self.log_bad_glue_source(-1, self.num_bg_lb_opt, self.num_bg_lb_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_igp_opt, self.num_bg_igp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_ogp_opt, self.num_bg_ogp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
                return new_bad_glues

        bad_glues_acg_fo,bad_glues_acg_g = self.check_algorithmic_conflicting_glues(temperature, len(new_bad_glues), num_bad_opt, threaded, weight_excess=ALGO_CONF_GLUES_WEIGHT_EXCESS)
        bad_glues_acg_fo = bad_glues_acg_fo*ALGO_CONF_GLUES_WEIGHT_EXTEND
        bad_glues_acg_g = bad_glues_acg_g*ALGO_CONF_GLUES_WEIGHT_EXTEND
        if not hasattr(self, 'num_bg_acg_orig_fo'): self.num_bg_acg_orig_fo = len(bad_glues_acg_fo)
        if not hasattr(self, 'num_bg_acg_opt_fo'): self.num_bg_acg_opt_fo = -1
        if not hasattr(self, 'num_bg_acg_orig_g'): self.num_bg_acg_orig_g = len(bad_glues_acg_g)
        if not hasattr(self, 'num_bg_acg_opt_g'): self.num_bg_acg_opt_g = -1
        self.num_bg_acg_most_recent_fo = len(bad_glues_acg_fo)
        self.num_bg_acg_most_recent_g = len(bad_glues_acg_g)
        self.log_bad_glue_source(len(bad_glues_acg_fo), self.num_bg_acg_opt_fo, self.num_bg_acg_orig_fo, '|')
        self.log_bad_glue_source(len(bad_glues_acg_g), self.num_bg_acg_opt_g, self.num_bg_acg_orig_g, '|')
        new_bad_glues.extend(bad_glues_acg_fo)
        new_bad_glues.extend(bad_glues_acg_g)

        if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt:
            if not CHECK_ONLY_ACG:
                self.log_bad_glue_source(-1, self.num_bg_lb_opt, self.num_bg_lb_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_igp_opt, self.num_bg_igp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_ogp_opt, self.num_bg_ogp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
            return new_bad_glues

        if not CHECK_ONLY_ACG:
            bad_glues_lb = self.check_lattice_binding(temperature, len(new_bad_glues), num_bad_opt, threaded)
            if not hasattr(self, 'num_bg_lb_orig'): self.num_bg_lb_orig = len(bad_glues_lb)
            if not hasattr(self, 'num_bg_lb_opt'): self.num_bg_lb_opt = -1
            self.num_bg_lb_most_recent = len(bad_glues_lb)
            self.log_bad_glue_source(len(bad_glues_lb), self.num_bg_lb_opt, self.num_bg_lb_orig, '|')
            new_bad_glues.extend(bad_glues_lb)

            if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt:
                self.log_bad_glue_source(-1, self.num_bg_igp_opt, self.num_bg_igp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_ogp_opt, self.num_bg_ogp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
                return new_bad_glues

            bad_glues_igp = self.check_input_pairs(temperature, len(new_bad_glues), num_bad_opt, threaded)
            if not hasattr(self, 'num_bg_igp_orig'): self.num_bg_igp_orig = len(bad_glues_igp)
            if not hasattr(self, 'num_bg_igp_opt'): self.num_bg_igp_opt = -1
            self.num_bg_igp_most_recent = len(bad_glues_igp)
            self.log_bad_glue_source(len(bad_glues_igp), self.num_bg_igp_opt, self.num_bg_igp_orig, '|')
            new_bad_glues.extend(bad_glues_igp)

            if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt:
                self.log_bad_glue_source(-1, self.num_bg_ogp_opt, self.num_bg_ogp_orig, '|')
                self.log_bad_glue_source(-1, self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
                return new_bad_glues

            bad_glues_ogp = self.check_output_pairs(temperature, len(new_bad_glues), num_bad_opt, threaded)
            if not hasattr(self, 'num_bg_ogp_orig'): self.num_bg_ogp_orig = len(bad_glues_ogp)
            if not hasattr(self, 'num_bg_ogp_opt'): self.num_bg_ogp_opt = -1
            self.num_bg_ogp_most_recent = len(bad_glues_ogp)
            self.log_bad_glue_source(len(bad_glues_ogp), self.num_bg_ogp_opt, self.num_bg_ogp_orig, '|')
            new_bad_glues.extend(bad_glues_ogp)

            if num_bad_opt > 0 and len(new_bad_glues) > num_bad_opt:
                self.log_bad_glue_source(-1, self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
                return new_bad_glues

            if check_tile_pairs:
                bad_glues_tp = self.check_tile_pairs_heuristic(temperature, len(new_bad_glues), num_bad_opt, threaded_tile_pairs, weight_excess=TILE_PAIRS_HEURISTIC_WEIGHT_EXCESS)
                bad_glues_tp = bad_glues_tp*TILE_PAIRS_HEURISTIC_WEIGHT_EXTEND
                if not hasattr(self, 'num_bg_tp_orig'): self.num_bg_tp_orig = len(bad_glues_tp)
                if not hasattr(self, 'num_bg_tp_opt'): self.num_bg_tp_opt = -1
                self.num_bg_tp_most_recent = len(bad_glues_tp)
                self.log_bad_glue_source(len(bad_glues_tp), self.num_bg_tp_opt, self.num_bg_tp_orig, '|')
                new_bad_glues.extend(bad_glues_tp)

        return new_bad_glues


    def update_cached_num_bad_glues(self):
        if hasattr(self, 'num_bg_acg_most_recent_fo'): self.num_bg_acg_opt_fo = self.num_bg_acg_most_recent_fo
        if hasattr(self, 'num_bg_acg_most_recent_g'): self.num_bg_acg_opt_g = self.num_bg_acg_most_recent_g
        if hasattr(self, 'num_bg_tss_most_recent'): self.num_bg_tss_opt = self.num_bg_tss_most_recent
        if hasattr(self, 'num_bg_igp_most_recent'): self.num_bg_igp_opt = self.num_bg_igp_most_recent
        if hasattr(self, 'num_bg_ogp_most_recent'): self.num_bg_ogp_opt = self.num_bg_ogp_most_recent
        if hasattr(self, 'num_bg_tp_most_recent'): self.num_bg_tp_opt = self.num_bg_tp_most_recent
        if hasattr(self, 'num_bg_lb_most_recent'): self.num_bg_lb_opt = self.num_bg_lb_most_recent

    def local_search_end_assign(self, temperature, out_filename, mutually_exclusive_tilename_pairs, seqs_start):
        '''Find assignment of unassigned glues to ends that avoids unwanted
        secondary structure in tiles, and also that maintains orthogonality
        between ends that could conflict (appear next to each other because
        they share an input glue on another input side), searching random
        permutations of ends.

        Uses stochastic local search heuristic.'''
        print('\n' + '*'*152 + '\nnumber unique glues:' + str(len(set([tile.glue(direction) for tile in self.tiles for direction in directions]))))

        print('threaded? ', THREADED)
        threaded = THREADED
        threaded_tile_pairs = THREADED_TILE_PAIRS

        self.num_bg_acg_opt_fo = -1
        self.num_bg_acg_opt_g = -1
        self.num_bg_tss_opt = -1
        self.num_bg_igp_opt = -1
        self.num_bg_ogp_opt = -1
        self.num_bg_tp_opt = -1
        self.num_bg_lb_opt = -1

        check_tile_pairs = (self.tiles[0].orth > 0 or self.tiles[0].orth_share > 0)
        if not check_tile_pairs:
            self.num_bg_tp_orig = -1

        self.assign_initial_ends(temperature, seqs_start)

        print('*'*152 + '\nformat of errors below indicating number of "bad glues" (glues causing problems meeting constraints) due to different sources: \n  <number in best glues found so far># <number in current iteration>/ <number in first iteration>\n' + '*'*152)

        print('Finding initial set of %s' % ('bad glues'))
        start = time.time()
        bad_glues_opt = self.find_bad_glues(temperature=temperature, num_bad_opt=-1, threaded=threaded, threaded_tile_pairs=threaded_tile_pairs, check_tile_pairs=check_tile_pairs)
        end = time.time()
        print('\n' + '*'*152 + ('\n** %.2f seconds for first call to find_bad_glue' % ((end - start))))

        self.update_cached_num_bad_glues()
        self.write_tile_sequences_to_file(out_filename, len(bad_glues_opt), 0)

        self.num_bad_glues_orig = len(bad_glues_opt)

        if check_tile_pairs:
            self.tile_pair_name2energy_opt = self.tile_pair_name2energy_most_recent

        iteration = 1
        print('*'*152 + '\nInitial number of bad glues %d:' % len(bad_glues_opt))
        picked_ends = list(set(glue.canonical_end() for glue in self.glues()))

        while len(bad_glues_opt) > 0:
            print('*'*152)
            iteration += 1

            glue = random.choice(bad_glues_opt)# DW: pick a bad glue
            ec = glue.end_constraint # DW: get EC for picked bad glue
            old_end = glue.canonical_end() # DW: get currently assigned seq for picked bad glue

            start=time.time()
            new_end_pos = select_new_end_pos(picked_ends, ec, temperature) # DW: randomly pick a new sequence
            end=time.time()
            # print '%.2f seconds to pick new end' % ((end - start))

            new_end = ec.ends[new_end_pos]

            # assign the new end
            glue.assign(new_end)
            self.tile_names_changed_from_opt = set(tile.name for tile in glue.tiles)
#             print 'glue changed: {}\ntile names changed from opt: {}'.format(glue, ', '.join(self.tile_names_changed_from_opt))

            start = time.time()
            new_bad_glues = self.find_bad_glues(temperature=temperature, num_bad_opt=len(bad_glues_opt),
                                                threaded=threaded, threaded_tile_pairs=threaded_tile_pairs,
                                                check_tile_pairs=check_tile_pairs)
            end = time.time()
            self.log_bad_glue_source(len(new_bad_glues), len(bad_glues_opt), self.num_bad_glues_orig, '')
            print('\niteration %d; %.2f seconds' % (iteration, (end - start)))

            num_worse = len(new_bad_glues) - len(bad_glues_opt)
#             print 'bad glues: %d in optimal sequences; %d if change made' % (len(bad_glues_opt), len(new_bad_glues))
#             print '%d bad glues previously' % len(bad_glues_opt)
#             print '%d bad glues if change made' % len(new_bad_glues)
            sys.stdout.flush()
            keep_change = (num_worse <= 0) # better glues (meaning, no worse): keep them!

            if not keep_change:
                glue.assign(old_end)
            else:
                # this is my hack to quickly swap out an old item for a new one
                ec.ends[new_end_pos] = old_end
                bad_glues_opt = new_bad_glues
                pos = picked_ends.index(old_end)
                picked_ends[pos] = new_end

                self.update_cached_num_bad_glues()

                # update energies of bad tile pairs
                # XXX: this must be done outside of check_tile_pairs_heuristic because we only want to update
                #      self.tile_pair_name2energy_opt if this is a new optimal sequence assignment
                if check_tile_pairs:
                    for ((t1name,t2name), energy) in self.tile_pair_name2energy_most_recent.items():
                        self.tile_pair_name2energy_opt[(t1name,t2name)] = energy
#                     old_energy = self.tile_pair_name2energy_opt[(t1name,t2name)]
#                     print '  updating tiles {:12s} and {:12s} from old energy {:2.2f} to new energy {:2.2f}'.format(t1name, t2name, old_energy, energy)

                self.write_tile_sequences_to_file(out_filename, len(bad_glues_opt), iteration)


    def write_tile_sequences_to_file(self, out_filename, num_bad_glues, iteration):
        ec = self.tiles[0].glue('N').end_constraint
        with open(out_filename, 'w') as f:
            print('** saving to %s with %d bad glues (weighted) **' % (out_filename, num_bad_glues))
            f.write('# saved at iteration %d\n' % iteration)
            #print 'num bad glues: {}'.format(len(bad_glues_opt))
            f.write('# %d bad glues\n' % num_bad_glues)
            f.write('# Algorithmic conflicting glues weight = %d' % ALGO_CONF_GLUES_WEIGHT_EXCESS)
            f.write(', tile secondary structure weight = %d' % TILE_SEC_STRUCT_WEIGHT_EXCESS)
            f.write(', tile pair heuristic weight = %d' % TILE_PAIRS_HEURISTIC_WEIGHT_EXCESS)
            f.write(', tile_sec_struct = ' + str(self.tiles[0].sec_struct))
            f.write(', tile_orth = ' + str(self.tiles[0].orth))
            f.write(', tile_orth_share = ' + str(self.tiles[0].orth_share))
            f.write(', orth_algorithmic_conflict = ' + str(ec.orth_algorithmic_conflict))
            f.write(', orth_algorithmic_conflict_generalized = ' + str(ec.orth_algorithmic_conflict_generalized))
            f.write(', lattice_binding_lower_threshold = ' + str(ec.lattice_binding_lower_threshold))
            f.write(', orth_colocated = ' + str(ec.orth_colocated) + '\n\n')
            max_len = max([len(tile.name) for tile in self.tiles])
            for tile in self.tiles:
                len_str = '%-' + ('%d' % max_len) + 's'
                f.write((('%s  ' % len_str) + '%s\n') % (tile, tile.sequence_spaced(include_biotins=True)))

    def assign_initial_ends(self, temperature, seqs_start):
        '''Assign ends to non-bound glues randomly.

        It will remove used ends from ec2unused_ends, as well as shuffle those lists.'''

        print('Assigning initial ends to glues {}'.format('based on sequences file' if seqs_start else 'randomly'))
        picked_ends = list(set(glue.canonical_end() for glue in self.glues_bound()))
        num_left = len(self.glues_unbound())
        for glue in self.glues_unbound():
#             print 'glues left to assign: %d,  assigning glue %s' % (num_left,glue)
            num_left -= 1
            if glue.assigned(): raise AssertionError("this shouldn't be reachable")
            ec = glue.end_constraint

#             start=time.time()
            if not seqs_start:
                new_end_pos = select_new_end_pos(picked_ends, ec, temperature)
            else:
                tilename = glue.canonical_tiles()[0].name
                direction = glue.canonical_direction()
                seq = seqs_start[tilename][direction]
                new_end_pos = ec.ends.index(seq)
#             end=time.time()
#             print 'time to find new end: %.2f' % (end-start)

            end = ec.ends[new_end_pos]
            glue.assign(end)
            last_pos = len(ec.ends) - 1
            ec.ends[new_end_pos], ec.ends[last_pos] = ec.ends[last_pos], ec.ends[new_end_pos]
            ec.ends.pop()
            picked_ends.append(end)


    def glues(self):
        '''Return all glues in all tiles in TileSet.'''
        return list(set(tile.glue(direction) for tile in self.tiles for direction in directions))

    def glues_unbound(self):
        '''Return all unbound glues in all tiles in TileSet.'''
        return [glue for glue in self.glues() if not glue.bound()]

    def glues_bound(self):
        '''Return all bound glues in all tiles in TileSet.'''
        return [glue for glue in self.glues() if glue.bound()]

    def __str__(self):
        return str(self.tiles)


def all_pairs_except(iterable, exclude):
    "Returns all pairs of elements (e1,e2) except those with exclude(e1,e2) true."
    def exclude_pack(args):
        e1=args[0]
        e2=args[1]
        return not exclude(e1,e2)
    return filter(exclude_pack, itertools.combinations_with_replacement(iterable, 2))

def eval_lattice_binding_energy(end1,end2,temperature):
    return min(sd.binding(end1+end2,wc(end2)+'TTTTT'+wc(end1),temperature),
               sd.binding(end1+end2,wc(end2)+'AAAAA'+wc(end1),temperature))

def eval_lattice_binding_energies(end_pairs, temperature, threaded):
    if threaded:
        results = [global_thread_pool.apply_async(eval_lattice_binding_energy, args=(end1, end2, temperature)) for (end1,end2) in end_pairs]
        energies = [result.get() for result in results]
    else:
        energies = [eval_lattice_binding_energy(end1, end2, temperature) for end1, end2 in end_pairs]
    return energies

def eval_colocated_end_pair(end1, end2, temperature):
#   return max(sd.hairpin(end1+end2, temperature), sd.hairpin(end1+'T'*4+end2, temperature))
    return max(sd.hairpin(end1+end2, temperature),
               min(sd.hairpin(end1+'T'*4+end2, temperature),
                   sd.hairpin(end1+'A'*4+end2, temperature)))

def eval_input_end_pair(end1, end2, temperature):
    return max(sd.hairpin(end1+end2, temperature))

def eval_colocated_end_pairs(end_pairs, temperature, threaded):
    if threaded:
        results = [global_thread_pool.apply_async(eval_colocated_end_pair, args=(end1, end2, temperature)) for (end1,end2) in end_pairs]
        energies = [result.get() for result in results]
    else:
        energies = [eval_colocated_end_pair(end1, end2, temperature) for end1, end2 in end_pairs]
    return energies

def share_complementary_glues(t1, t2):
    for direction in directions:
        if t1.glue(direction) == t2.glue(opposite(direction)):
            return True
    return False

LOG_BAD_END = False
def log_bad_end(reason):
    if LOG_BAD_END:
        sys.stdout.write(reason)
        sys.stdout.flush()

def select_new_end_pos(picked_ends, ec, temperature, threaded=True):
    num_searched = 0
    if LOG_BAD_END: sys.stdout.write('.')
    if LOG_BAD_END: sys.stdout.flush()

    picked_ends_len = dict()
    picked_ends_len[10] = [end for end in picked_ends if len(end) == 10]
    picked_ends_len[11] = [end for end in picked_ends if len(end) == 11]

    for end_pos in random_iter(range(len(ec.ends))):
        end = ec.ends[end_pos]
        num_searched += 1
        if wc(end) in picked_ends:
            log_bad_end('wc_')
            continue
        if end in picked_ends:
            log_bad_end('used_')
            continue
        if (ec.lowPF is not None) and ec.lowPF > 0 and (ec.highPF is not None) and ec.highPF > 0 and not sd.domain_equal_strength(end,temperature,ec.lowPF,ec.highPF):
            log_bad_end('eq_')
            continue
        if ec.domain_indv_sec_struct > 0 and not sd.domain_no_sec_struct(end,temperature,ec.domain_indv_sec_struct,threaded):
            log_bad_end('idv%.1f,%.1f_' % (sd.hairpin(end,temperature), sd.hairpin(sd.wc(end),temperature)))
            continue
        if ec.domain_pair_sec_struct > 0:
            ave = ec.domain_pair_sec_struct_ave if ec.domain_pair_sec_struct_ave > 0 else -1
            good = sd.domain_pairwise_concatenated_no_sec_struct(end,picked_ends,temperature,ec.domain_pair_sec_struct,ave,threaded)
            if not good:
                log_bad_end('pair_')
                continue
        if ec.orth_any > 0:
            ave = ec.orth_any_ave if ec.orth_any_ave > 0 else -1
            good = sd.domain_orthogonal(end,picked_ends,temperature,ec.orth_any,ave,threaded)
            if not good:
                log_bad_end('orth_')
                continue
        if ec.hamming > 0:
            good = _eval_pairs_hamming(end, picked_ends_len[len(end)], ec.hamming)
            if not good:
                log_bad_end('ham_')
                continue
        if LOG_BAD_END: sys.stdout.write('.\n')
        if LOG_BAD_END: sys.stdout.flush()
        return end_pos
    raise ValueError('no more sequences to search')
    return None

def hamming(s1,s2):
    '''Return minimum Hamming distance between s1 and s2, against all shifts if one is longer.'''
    if len(s1) == len(s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    if len(s1) > len(s2): s1,s2 = s2,s1
    excess = len(s2) - len(s1)
    return min( sum(c1 != c2 for c1, c2 in zip(s1, s2[start:start+len(s1)]))
                for start in range(excess+1))

def hamming_with_wc(end1,end2):
    return min(hamming(end1,end2), hamming(end1,wc(end2)))

def _eval_pairs_hamming(end, picked_ends, hamming_thres):
    return min(hamming_with_wc(end, picked_end) for picked_end in picked_ends) >= hamming_thres

def _eval_pairs(endpairs, temperature, parallel):
    if parallel:
        results = [global_thread_pool.apply_async(sd.binding, args=(end1, end2, temperature)) for (end1,end2) in endpairs]
        energies = [result.get() for result in results]
    else:
        energies = []
        for end1,end2 in endpairs:
            energies.append(sd.binding(end1, end2, temperature))
    return energies

def eval_end_pair(end1,end2,temperature,parallel=True):
    end1wc = wc(end1)
    end2wc = wc(end2)
    if end1 == end2:
        energy1_2, energy1wc_2wc = _eval_pairs([(end1,end2), (end1wc, end2wc)], temperature, parallel)
        energy1wc_2 = energy1_2wc = 0.0
    elif end1 == end2wc:
        energy1wc_2, energy1_2wc = _eval_pairs([(end1wc,end2), (end1, end2wc)], temperature, parallel)
        energy1_2 = energy1wc_2wc = 0.0
    else:
        energy1_2, energy1wc_2wc, energy1wc_2, energy1_2wc = _eval_pairs([(end1,end2), (end1wc,end2wc), (end1wc,end2), (end1, end2wc)], temperature, parallel)
    return energy1_2, energy1_2wc, energy1wc_2, energy1wc_2wc

def prefilter_ends(end_constraints, temperature, three_letter_code, three_letter_code_exceptions, user_assigned_glue_sequences):
    '''Filter ends that can be quickly filtered using NumPy.

    Place a member variable called ends into each EndConstraint, a list of
    potential ends satisfying the constraints that can be quickly checked.

    If three_letter_code is True, then parity 1 tiles will be from the alphabet
    ('A','C','T') and parity 0 tiles will be from the alphabet ('A','G','T')'''
    forbidden_subs = ['%s%s%s%s' % (a, b, c, d) for a in ('G', 'C')
                                                for b in ('G', 'C')
                                                for c in ('G', 'C')
                                                for d in ('G', 'C')]
    all_lengths = { ec.length for ec in end_constraints }
    ends_of_length = { length:dsd.DNASeqList(length=length) for length in all_lengths }

    for (length, ends) in ends_of_length.items():
#         ec = [ec for ec in end_constraints if ec.length == length][0]
#         ends_of_length[length] = ends
        print('Removing length-%d domains containing {C,G}^4' % length)
        ends = ends.filter_substring(forbidden_subs)
        ends_of_length[length] = ends
#         ends_lst = ends.toList()
#         for end in ends_lst:
#             if 'GGGG' in end:
#                 print 'ALERT %s' % end
#             if 'CCCC' in end:
#                 print 'ALERT %s' % end

    uags = {} # dictionary of user assigned glue sequences
    #print user_assigned_glue_sequences
    for length in all_lengths:
        uags[length] = [ seq for (_,_,seq) in user_assigned_glue_sequences if len(seq)==length] + \
                [ sd.wc(seq) for (_,_,seq) in user_assigned_glue_sequences if len(seq)==length]

    for ec in end_constraints:
        ends = ends_of_length[ec.length]
        if ec.endGC:
            print('Removing length-%d domains that end in A or T' % length)
            ends = ends.filter_endGC()
        if ec.endAT:
            print('Removing length-%d domains that end in G or C' % length)
            ends = ends.filter_endAT(gc_near_end=False)
        if three_letter_code:
            base_to_remove = 'G' if ec.glues[0].input_parity() else 'C'
            print("Removing length-%d domains that contain > %d %s's" % (ec.length,three_letter_code_exceptions,base_to_remove))
            ends = ends.filter_base_count(base_to_remove, 0, three_letter_code_exceptions)
        if ec.biotin_direction:
            biotin_direction_parity = ec.glues[0].get_parity(ec.biotin_direction)
            biotin_pos = biotin_end_pos(ec.biotin_direction, biotin_direction_parity)
            if is_canonical(ec.biotin_direction):
                base = 'T'
                biotin_canonical_pos = biotin_pos
                canonical_direction = ec.biotin_direction
            else:
                base = 'A'
                biotin_canonical_pos = ec.length - 1 - biotin_pos
                canonical_direction = opposite(ec.biotin_direction)
            print(('Removing length-%d domains lacking %s at pos %d on %s-facing domain (to place biotin at position %d on %s-facing domain)'
                  % (ec.length,base, biotin_canonical_pos, canonical_direction, biotin_pos, ec.biotin_direction)))
            ends = ends.filter_base_at_pos(biotin_canonical_pos, base)
#         print 'before filtering based on energy: %d' % ends.numseqs
        ends = ends.filter_energy(ec.lowDG, ec.highDG, temperature)
        print(('NN binding energy: %d length-%d seqs in (%.1f,%.1f) %s'
               % (ends.numseqs, ec.length, ec.lowDG, ec.highDG,
                  '' if not ec.biotin_direction else '[biotin pos ' + str(biotin_pos) + ']')))
        ec.ends = ends.toList()

        sys.stdout.write('Removing from the set of length-{} domains, if they are present, the following domains that are already assigned to glues by user: '.format(ec.length) + ', '.join(s for s in uags[ec.length]) +'. ')
        #print str(len(ec.ends))+ ' number ends before'
        tmp_ends_before = ec.ends
        ec.ends = [end for end in ec.ends if end not in uags[ec.length]]
        if len(ec.ends) == 0:
            raise ValueError('error: no ends of length %d found matching requested parameters' % ec.length)
        print(str(len(ec.ends))+ ' ends after removal of: ' + str(', '.join(list(set(tmp_ends_before).difference(set(ec.ends))))
           + ' None'*(', '.join(list(set(tmp_ends_before).difference(set(ec.ends))))=='') ))

try:
    _zip_longest = itertools.zip_longest # python3
except:
    _zip_longest = itertools.izip_longest # python2

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return _zip_longest(fillvalue=fillvalue, *args)

def alphabet_to_use(three_letter_code, parity, direction):
    '''Return tuple of alphabet to be used for glue in given direction on tile of
    given parity.

    Note that this refers to the alphabet used for the CANONICAL direction, which
    may be the opposite of direction.'''
    if not parity in (0,1):
        raise ValueError('parity must be 0 or 1, cannot be %s' % parity)
    if not direction in directions:
        raise ValueError('direction must be in %s, cannot be %s' % (directions, direction))
    if not three_letter_code:
        return ('A','C','G','T')
    if (parity == 1 and is_canonical(direction)) or (parity == 0 and not is_canonical(direction)):
        return ('A','C','T')
    else:
        return ('A','G','T')

def read_tiles_from_file(filename):
    '''Example of format:

    # "global" energy interval for glues (not set by glue_strength)
    lowDG = 9.8
    highDG = 10.2

    # temperature at which to evaluate energies
    temperature = 53.0

    # threshold for allowed energy of secondary structure
    tile_sec_struct = 2.0

    tile_pair_sec_struct = 7.0

    # threshold for allowed energy of interaction between orthogonal glues
    # that might be co-located during a strength-1 binding event
    orth_algorithmic_conflict = 2.0

    # threshold for allowed energy of interaction between orthogonal glues
    # that will not be co-located during a strength-1 binding event
    orth_any = 4.0

    # orth_any is worst-case; all pairs of sequences must be below
    # orth_any_ave is average orthogonality of proposed new domain with all
    # currently used domains
    orth_any_ave = 3.0

    # this is the limit on individual secondary structure of a single domain
    domain_indv_sec_struct = 0.8

    # this is the limit on secondary structure of any two domains concatenated
    domain_pair_sec_struct = 1.2

    # if three_letter_code == true, then parity 1 tiles will use a 3-letter code
    #  A,C,T, and parity 0 tiles will use A,G,T
    three_letter_code = False


    endGC = True

    # XOR tiles
    tiles = [ { 'name':'xor00_0', 'N':"0",  'W':"0",  'S':"0'", 'E':"0'", 'parity':0 },
              { 'name':'xor01_0', 'N':"0",  'W':"1",  'S':"1'", 'E':"1'", 'parity':0 },
              { 'name':'xor10_0', 'N':"1",  'W':"0",  'S':"1'", 'E':"1'", 'parity':0 },
              { 'name':'xor11_0', 'N':"1",  'W':"1",  'S':"0'", 'E':"0'", 'parity':0 },
              { 'name':'xor00_1', 'N':"0'", 'W':"0'", 'S':"0",  'E':"0",  'parity':1 },
              { 'name':'xor01_1', 'N':"0'", 'W':"1'", 'S':"1",  'E':"1",  'parity':1 },
              { 'name':'xor10_1', 'N':"1'", 'W':"0'", 'S':"1",  'E':"1",  'parity':1 },
              { 'name':'xor11_1', 'N':"1'", 'W':"1'", 'S':"0",  'E':"0",  'parity':1 } ]

    # glues that need different strengths from lowDG and highDG
    glue_strength_constraints = [
                                 ( 10.5,11.5, [ ("0'","NS") ] ),
                                 (  9.5, 9.8, [ ("0", "NS") ] )
                                ]
    # glues that need biotins
    # (hence must put an internal biotinylated dT in end)
    glue_biotins = [ ("1","S") ]

    # sequences to bind to glues
    glue_sequences = [ ("0",  "W", "TTGAGGAGAG"),
                       ("0'", "W", "GTGTAGTAGGC")]

    '''
    import imp
    mod = imp.load_source('',filename)

    lowDG = mod.lowDG
    highDG = mod.highDG
    temperature = mod.temperature
    tile_sec_struct = mod.tile_sec_struct
    tile_orth = mod.tile_orth
    tile_orth_share = mod.tile_orth_share
    orth_algorithmic_conflict = mod.orth_algorithmic_conflict
    orth_colocated = mod.orth_colocated
    orth_any = mod.orth_any
    orth_any_ave = mod.orth_any_ave
    three_letter_code = mod.three_letter_code
    three_letter_code_exceptions = mod.three_letter_code_exceptions
    endGC = mod.endGC
    endAT = mod.endAT
    domain_indv_sec_struct = mod.domain_indv_sec_struct
    domain_pair_sec_struct = mod.domain_pair_sec_struct
    domain_pair_sec_struct_ave = mod.domain_pair_sec_struct_ave
    hamming = mod.hamming
    orth_algorithmic_conflict_generalized = mod.orth_algorithmic_conflict_generalized
    lattice_binding_lower_threshold = mod.lattice_binding_lower_threshold
    user_assigned_glue_sequences = mod.glue_sequences

    detect_nondeterminism = mod.detect_nondeterminism if hasattr(mod, "detect_nondeterminism") else True

    glue_strength = collections.defaultdict(lambda: (lowDG, highDG))
    if hasattr(mod, 'glue_strength_constraints'):
        for (lowDG_e,highDG_e,glues) in mod.glue_strength_constraints:
            for label,axis in glues:
                if axis not in ['NS','EW']:
                    raise ValueError('axis is {}; should be one of "NS" or "EW"'.format(axis))
                glue_strength[(label,axis)] = (lowDG_e,highDG_e)

    glue_biotin = collections.defaultdict(lambda: None)
    if hasattr(mod, 'glue_biotins'):
        for (label,direction) in mod.glue_biotins:
            if direction not in directions: raise ValueError('"%s" is not a valid direction' % direction)
            axis = direction2axis(direction)
            if (label,axis) in glue_biotin:
                raise ValueError('error: glues %s on axis %s already assigned a biotin direction of %s \n cannot assign multiple directions'
                                 % (label, axis, glue_biotin[(label,axis)]))
            else:
                glue_biotin[(label,axis)] = direction


    dna_alphabet_pattern = re.compile('(A|C|G|T)+')

    glue_seq = dict()
    if hasattr(mod, 'glue_sequences'):
        for (label,direction,seq) in mod.glue_sequences:
            axis = direction2axis(direction)
            if not is_canonical(direction):
                seq = wc(seq)
            if not dna_alphabet_pattern.match(seq):
                raise ValueError('{} is not a DNA string'.format(seq))
            glue_seq[(label, axis)] = seq

    if hasattr(mod, 'mutually_exclusive_tilename_pairs'):
        mutually_exclusive_tilename_pairs = mod.mutually_exclusive_tilename_pairs
    else:
        mutually_exclusive_tilename_pairs = []

    if not hasattr(mod, 'tiles'):
        raise ValueError('''must define a list of dicts called tiles, e.g.,

  [ { 'name':'xor00_0', 'N':"0",  'W':"0",  'S':"0'", 'E':"0'", 'parity':0 },
    { 'name':'xor01_0', 'N':"0",  'W':"1",  'S':"1'", 'E':"1'", 'parity':0 } ]''')


    tile_dicts = mod.tiles
    tiles = []

    print('processing read tile descriptions')

    for tile_dict in tile_dicts:
        n_label = tile_dict['N']
        s_label = tile_dict['S']
        e_label = tile_dict['E']
        w_label = tile_dict['W']
        name = tile_dict['name']
        parity = tile_dict['parity']

        glues = dict()
        for (label,direction) in zip([n_label,s_label,e_label,w_label] , ['N','S','E','W']):
            axis = direction2axis(direction)
            lowDG_e,highDG_e = glue_strength[(label,axis)]
            biotin_direction = glue_biotin[(label,axis)]
            length = len_end(direction, parity)
            ec = EndConstraint.factory(lowDG=lowDG_e, highDG=highDG_e,
                length=length, biotin_direction=biotin_direction, endGC=endGC, endAT=endAT,
                orth_any=orth_any, orth_any_ave=orth_any_ave,
                orth_algorithmic_conflict=orth_algorithmic_conflict, orth_colocated=orth_colocated,
                domain_indv_sec_struct=domain_indv_sec_struct,
                domain_pair_sec_struct=domain_pair_sec_struct, domain_pair_sec_struct_ave=domain_pair_sec_struct_ave,
                tile_sec_struct=tile_sec_struct,
                hamming = hamming,
                orth_algorithmic_conflict_generalized = orth_algorithmic_conflict_generalized,
                lattice_binding_lower_threshold=lattice_binding_lower_threshold)
            glue = glues[direction] = Glue.factory(label=label,axis=axis,end_constraint=ec)
            if (label, axis) in glue_seq:
                canonical_seq = glue_seq[(label, axis)]
                if not glue.bound():
                    glue.bind(canonical_seq)
                else:
                    if canonical_seq != glue.canonical_end():
                        print('Cannot bind sequence %s to glue %s because it is already bound to sequence %s' % (canonical_seq, glue, glue.canonical_end()))
                end = glue.get_end(direction, include_biotins=False)
                if direction == biotin_direction:
                    biotin_pos = biotin_end_pos(direction, parity)
                    if end[biotin_pos] != 'T':
                        raise ValueError('to include internal biotin on glue %s in direction %s, base at position %d must be T, but instead is %s' %
                                 (glue, direction, biotin_pos, end[biotin_pos]))
        if 'tile_sec_struct' in tile_dict:
            tile_sec_struct_local = tile_dict['tile_sec_struct']
        else:
            tile_sec_struct_local = tile_sec_struct
        if 'tile_orth' in tile_dict:
            tile_orth_local = tile_dict['tile_orth']
        else:
            tile_orth_local = tile_orth
        if 'frna_share' in tile_dict:
            tile_orth_share_local = tile_dict['tile_orth_share']
        else:
            tile_orth_share_local = tile_orth_share
        tile = Tile(name=name, n=glues['N'], s=glues['S'], e=glues['E'], w=glues['W'], parity=parity, sec_struct=tile_sec_struct_local, orth=tile_orth_local, orth_share=tile_orth_share_local)
        tiles.append(tile)

    for tile in tiles:
        for glue in [tile.glue(direction) for direction in directions]:
            if tile not in glue.tiles:
                raise ValueError('ERROR: tile {} not in glue.tiles for glue={}'.format(tile.name, glue.label))

    print('done processing tiles')

    tileset = TileSet(tiles, detect_nondeterminism=detect_nondeterminism, mutually_exclusive_tilename_pairs=mutually_exclusive_tilename_pairs)

    print('done creating tile set')

    if hasattr(mod, 'glue_biotins'):
        for (label,direction) in mod.glue_biotins:
            if label not in [glue.label for glue in tileset.glues()]:
                raise ValueError('%s is not a valid glue label' % label)

    if hasattr(mod, 'glue_strength_constraints'):
        for (lowDG,highDG,labels) in mod.glue_strength_constraints:
            for (label,axis) in labels:
                if label not in [glue.label for glue in tileset.glues()]:
                    raise ValueError('%s is not a valid glue label' % label)

    if hasattr(mod, 'glue_sequences'):
        for (label,direction,seq) in mod.glue_sequences:
            if label not in [glue.label for glue in tileset.glues()]:
                raise ValueError('%s is not a valid glue label' % label)

    tile_names = [tile.name for tile in tiles]
    if hasattr(mod, 'mutually_exclusive_tilename_pairs'):
        for t1,t2 in mod.mutually_exclusive_tilename_pairs:
            if t1 not in tile_names:
                raise ValueError('%s is not a valid tile name' % t1)
            if t2 not in tile_names:
                raise ValueError('%s is not a valid tile name' % t2)

    return tileset, temperature, three_letter_code, three_letter_code_exceptions, mutually_exclusive_tilename_pairs, user_assigned_glue_sequences

def read_starting_seqs(seqs_filename):
    seqs_start = dict()
    with open(seqs_filename, 'r') as f:
        lines = f.readlines()
        lines = [line[:line.find('#')].strip() for line in lines]
        lines = [line for line in lines if len(line) > 0]
        for line in lines:
            name,s,w,n,e = line.split()
            s = s.replace('/iBiodT/', 'T')
            w = w.replace('/iBiodT/', 'T')
            n = n.replace('/iBiodT/', 'T')
            e = e.replace('/iBiodT/', 'T')
            seqs_start[name] = {'S':s, 'W':w, 'N':n, 'E':e}

    return seqs_start

def check_for_nupack():
    print('Checking for NUPACK... ', end=' ')
    try:
        sd.duplex("ACGT", 53)
    except:
        print('NUPACK is not installed correctly. Please install it and ensure that pfunc can be called from the command line, and that NUPACKHOME is appropriately set.')
        sys.exit(-1)
    print('NUPACK is installed correctly')


def check_for_viennarna():
    print('Checking for ViennaRNA... ', end=' ')
    try:
        sd.RNAduplex_multiple([("ACGT","TGCA")], 53)
    except:
        print('Vienna RNA is not installed correctly. Please install it and ensure that RNAduplex can be called from the command line.')
        sys.exit(-1)
    print('Vienna RNA is installed correctly')

import argparse

if __name__ == "__main__":
    check_for_viennarna()
    check_for_nupack()
    parser = argparse.ArgumentParser(description='Design SST strands to implement an abstract aTAM tile set')
    parser.add_argument('-p','--params', help='name of input parameter file', required=True)
    parser.add_argument('-o','--out', help='name of output file', required=True)
    parser.add_argument('-s','--seqs', help='name of file (in format of output file) indicating which sequences to start with', required=False)
    args = vars(parser.parse_args())

    in_filename = args['params']
    out_filename = args['out']
    seqs_filename = args.get('seqs')

    seqs_start = read_starting_seqs(seqs_filename) if seqs_filename else None

    print('number of processes in global thread pool: %d' % global_thread_pool._processes)

    print('Reading tiles from %s' % in_filename)

    tileset, temperature, tlc, tlce, mutually_exclusive_tilename_pairs, usgs = read_tiles_from_file(in_filename)

    prefilter_ends(tileset.end_constraints_unbound(), temperature, tlc, tlce, usgs)

    tileset.local_search_end_assign(temperature=temperature,
                                    out_filename=out_filename,
                                    mutually_exclusive_tilename_pairs=mutually_exclusive_tilename_pairs,
                                    seqs_start=seqs_start)

    print()
    for tile in tileset.tiles:
        print('%5s %s' % (tile, tile.sequence_spaced(include_biotins=True)))
