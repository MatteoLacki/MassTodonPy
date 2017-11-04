try:
    import cPickle as pickle
except ImportError:
    import pickle
import re
from linearCounter.linearCounter import linearCounter as lCnt
from collections import defaultdict


class NegativeAtomCount(Exception):
    pass


def atomCnt2string(atomCnt):
    """Translate a dictionary of atom counts into a uniquely defined string.

    Parameters
    ----------
    atomCnt : Counter
        The chemical formula counter.

    Returns
    -------
    out : str
        A chemical formula string.
    """
    keys = atomCnt.keys()
    return "".join(el+str(atomCnt[el]) for el in sorted(keys))


def standardize_modifications(modifications):
    """Standardize modifications so that they meet
    the internal nomenclature scheme.

    Parameters
    ----------
    atomCnt : Counter
        The chemical formula counter.

    Returns
    -------
    out : defaultdict
        The atomic modifications.

    Notes
    -----
    It was easier for me to think of an amino acid as if it was composed
    out of three bricks: the left one, the center one, and the right one.
    The left one corresponds to the group with nitrogen,
    the center one - to the alpha carbon (including the side chain),
    and the right one - to the other carbon atom.
    """
    backboneAtom2aaNomen = {'N': 'L', 'Calpha': 'C', 'C': 'R'}
    R = defaultdict(lambda: defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        match = re.match(r"([a-z]+)([0-9]+)", tag, re.I)
        if match:
            aa, aa_idx = match.groups()
            aa_idx = int(aa_idx) - 1
        R[aa_idx][backboneAtom2aaNomen[aa]] = lCnt(atomCnt)
    return R


def prolineBlockedFragments(fasta):
    """Checks for c-z fragments that cannot occur because of proline.

    Parameters
    ----------
    fasta : str
        The fasta of the studied molecular species.

    Returns
    -------
    blocked : set
        A set of fragment names that cannot occur.
        Always contains the 'c0', which is too small to observe.
    """
    blocked = set('c0')
    for i, f in enumerate(fasta):
        if f == 'P':
            blocked.add('c' + str(i))
            blocked.add('z' + str(len(fasta)-i))
    return blocked


class Formulator(object):
    def __init__(self, fasta, Q, data_path,
                 distance_charges=5,
                 modifications=lCnt()):
        self.Q = Q
        self.fasta = fasta
        self.d_charges = distance_charges
        self.modifications = standardize_modifications(modifications)
        with open(data_path+"amino_acids.txt", "rb") as f:
            self.bricks = pickle.load(f)

    def getBrick(self, aaPart, aa, aaNo):
        brick = self.bricks[aa][aaPart] + self.modifications[aaNo][aaPart]
        if any(brick[elem] < 0 for elem in brick):
            raise NegativeAtomCount("Attention: your modification had an unexpected effect.\
            Part of your molecule now has negative atom count.\
            Bear that in mind while publishing your results.")
        return brick

    def protonate(self, frag):
        a, b, c = {'p': (1, 0, 1),
                   'c': (0, -1, 0),
                   'z': (0, 0, 1)}[frag]

        for q in range(1, self.Q+a):
            for g in range(b, self.Q-q+c):
                yield (q, g)

    def make_fragments(self):
        superAtoms = []
        sA = lCnt()
        for aaNo, aa in enumerate(self.fasta):
            sA += self.getBrick('L', aa, aaNo)
            superAtoms.append(sA)
            sA = self.getBrick('C', aa, aaNo) + self.getBrick('R', aa, aaNo)
        sA += lCnt({'O': 1, 'H': 1})
        superAtoms.append(sA)
        superAtoms[0] += lCnt({'H': 1})

        N = len(superAtoms)
        N_fasta = len(self.fasta)

        # precursor molecule
        precursor = ('precursor', atomCnt2string(sum(superAtoms)), N_fasta)
        unprotonated_formulas = [precursor]

        # exclude proline fragmented daughter ions
        blockedFragments = prolineBlockedFragments(self.fasta)

        # c fragments
        # Adding one extra hydrogen to meet the definition of a c fragment.
        cFrag = lCnt({'H': 1})
        for i in range(N-1):
            cFrag += superAtoms[i]
            cFrag_tmp = lCnt(cFrag)
            frag_type = 'c'+str(i)
            if frag_type not in blockedFragments and not i == 0:
                c_fragment = (frag_type, atomCnt2string(cFrag_tmp), i)
                unprotonated_formulas.append(c_fragment)

        # z fragments
        zFrag = lCnt()
        for i in range(1, N):
            zFrag += superAtoms[N-i]
            zFrag_tmp = lCnt(zFrag)
            frag_type = 'z'+str(i)
            if frag_type not in blockedFragments:
                z_fragment = (frag_type, atomCnt2string(zFrag_tmp), i)
                unprotonated_formulas.append(z_fragment)

        # outcome variable
        formulas = []

        for molType, atomCnt_str, sideChainsNo in unprotonated_formulas:
            for q, g in self.protonate(molType[0]):
                potentialChargesNo = sideChainsNo // self.d_charges
                if sideChainsNo % self.d_charges > 0:
                    potentialChargesNo += 1
                    # +0000 +0000 00+  at most 3 charges
                if potentialChargesNo >= q:
                    formulas.append((molType, atomCnt_str, sideChainsNo, q, g))

        return formulas


import json
import numpy as np
data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"

with open(data_path + "substanceP.json","rb") as f:
    mol = json.load(f)
mol["spectrum"] = tuple( np.array(d) for d in mol["spectrum"] )


F = Formulator(fasta=mol["fasta"],
               Q=mol["Q"],
               data_path=data_path,
               modifications=mol["modifications"])

formulas = F.make_fragments()

with open(data_path+"amino_acids.txt", "rb") as f:
    bricks = pickle.load(f)
