# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

from linearCounter import linearCounter as lCnt
from itertools import chain
from protonations import protonate
import  pkg_resources
from collections import defaultdict
import re
try:
   import cPickle as pickle
except:
   import pickle


data_path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
bricks  = pickle.load(open(data_path+'amino_acids.txt', 'rb'))

def countIsNegative(atomCnt):
    '''Check if any element of a dictionary is a negative number.'''
    return any( atomCnt[elem]<0 for elem in atomCnt )


def atomCnt2string(atomCnt):
    '''Translate a dictionary of atom counts into a uniquely defined string.'''
    keys = atomCnt.keys()
    keys.sort()
    return "".join( el+str(atomCnt[el]) for el in keys )


def standardize(modifications):
    '''Standardize modifications so that they meet the internal nomenclature scheme.'''
    backboneAtom2aaNomen = {'N':'L', 'Calpha':'C', 'C':'R'}
    R = defaultdict(lambda:defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        match = re.match(r"([a-z]+)([0-9]+)", tag, re.I)
        if match:
            aa, aa_idx = match.groups()
            aa_idx = int(aa_idx) - 1
        R[aa_idx][ backboneAtom2aaNomen[aa] ] = lCnt(atomCnt)
    return R


def prolineBlockedFragments(fasta):
    '''Checks which c-z fragments cannot occur.'''
    blocked = set('c0')
    for i, f in enumerate(fasta):
        if f=='P':
            blocked.add( 'c' + str(i) )
            blocked.add( 'z' + str( len(fasta)-i ) )
    return blocked

def make_cz_fragments(fasta, modifications):
    '''Prepares the precursor and the c and z fragments atom counts.'''
    def getBrick(aaPart, aa):
        brick = bricks[aa][aaPart] + modifications[aaNo][aaPart]
        if countIsNegative(brick):
            print("Attention: your modification has an unexpected effect. Part of your molecule now has negative atom count. Bear that in mind while publishing your results.")
        return brick

    superAtoms = []
    sA = lCnt()
    for aaNo, aa in enumerate(fasta):
        sA += getBrick('L', aa)
        superAtoms.append( sA )
        sA = getBrick('C', aa) + getBrick('R', aa)
    sA += lCnt({'O':1,'H':1})
    superAtoms.append(sA)
    superAtoms[0] += lCnt({'H':1})
    N = len(superAtoms)

    def getPrecursor():
        precursor = sum(superAtoms)
        yield ('precursor', atomCnt2string(precursor), len(fasta) )

    blockedFragments = prolineBlockedFragments(fasta)

    def getCfrags():
        cFrag = lCnt({'H':1}) # Adding one extra hydrogen to meet the definition of a c fragment.
        for i in range(N-1):
            cFrag += superAtoms[i]
            cFrag_tmp = lCnt(cFrag)
            frag_type = 'c'+str(i)
            if not frag_type in blockedFragments and not i == 0:
                yield (frag_type, atomCnt2string(cFrag_tmp), i)

    def getZfrags():
        zFrag = lCnt()
        for i in range(1,N):
            zFrag += superAtoms[N-i]
            zFrag_tmp = lCnt(zFrag)
            frag_type = 'z'+str(i)
            if not frag_type in blockedFragments:
                yield (frag_type, atomCnt2string(zFrag_tmp), i)

    return getPrecursor, getCfrags, getZfrags
#TODO It seems very strange to return these functions. Inspect it later on.

class Formulator(object):
    def __init__(self, fasta, Q, modifications={} ):
        self.fasta  = fasta
        self.Q      = Q
        self.modifications = modifications

class CZformulator(Formulator):
    def __init__(self, fasta, Q, modifications={} ):
        super(CZformulator,self).__init__(fasta, Q, modifications)
        self.precs, self.cfrags, self.zfrags = make_cz_fragments(fasta, modifications)

                            # aaPerOneCharge
    def makeMolecules(self, distanceBetweenCharges=5):
        '''Generate possible molecules in c/z fragmentation.

        Returns: tuples ( type_of_mol, mol_formula, amino_acids_no, charge, quenched_charge ).'''
        for molType, atomCnt_str, sideChainsNo in chain( self.precs(), self.cfrags(), self.zfrags() ):
            for q,g in protonate( self.Q, molType[0] ):
                potentialChargesNo = sideChainsNo / distanceBetweenCharges
                if sideChainsNo % distanceBetweenCharges > 0:
                    potentialChargesNo += 1 # +0000 +0000 00+  at most 3 charges
                if potentialChargesNo >= q:
                    yield molType, atomCnt_str, sideChainsNo, q, g

class CZformulator_qg_competition(CZformulator):
    def makeMolecules(self, distanceBetweenCharges=5):
        '''Generate possible molecules in c/z fragmentation. Take into account that q and g charges compete for the bloody places.

        Returns: tuples ( type_of_mol, mol_formula, amino_acids_no, charge, quenched_charge ).'''
        for molType, atomCnt_str, sideChainsNo in chain( self.precs(), self.cfrags(), self.zfrags() ):
            for q,g in protonate( self.Q, molType[0] ):
                if g >= 0:
                    totalCharges = q+g
                else:
                    totalCharges = q
                potentialChargesNo = sideChainsNo / distanceBetweenCharges
                if sideChainsNo % distanceBetweenCharges > 0:
                    potentialChargesNo += 1 # +0000 +0000 00+  at most 3 charges
                if potentialChargesNo >= totalCharges:
                    yield molType, atomCnt_str, sideChainsNo, q, g


def make_formulas(
        fasta,
        Q,
        frag_type ='cz',
        modifications={}
    ):
    '''Generate all possible fragments given a Roepstorf Scheme [or its generalization].
    '''
    modifications   = standardize(modifications)
    formClass       = {
        'cz':CZformulator,
        'cz_qg_competition':CZformulator_qg_competition
    }[frag_type](fasta, Q, modifications )
    return formClass
