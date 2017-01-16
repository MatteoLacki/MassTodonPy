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
from bricks import makeBricks
from misc import standardize, countIsNegative, atomCnt2string

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
    bricks = makeBricks()
    def getBrick(aaPart):
        brick = bricks[aa][aaPart] + modifications[aaNo][aaPart]
        if countIsNegative(brick):
            print("Attention: your modification has an unexpected effect. Part of your molecule now has negative atom count. Bear that in mind while publishing your results.")
        return brick

    superAtoms = []
    sA = lCnt()
    for aaNo, aa in enumerate(fasta):
        sA += getBrick('L')
        superAtoms.append( sA )
        sA = getBrick('C') + getBrick('R')
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
            fragType = 'c'+str(i)
            if not fragType in blockedFragments and not i == 0:
                yield (fragType, atomCnt2string(cFrag_tmp), i)

    def getZfrags():
        zFrag = lCnt()
        for i in range(1,N):
            zFrag += superAtoms[N-i]
            zFrag_tmp = lCnt(zFrag)
            fragType = 'z'+str(i)
            if not fragType in blockedFragments:
                yield (fragType, atomCnt2string(zFrag_tmp), i)

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

    def makeMolecules(self, aaPerOneCharge=5):
        for molType, atomCnt_str, sideChainsNo in chain( self.precs(), self.cfrags(), self.zfrags() ):
            for q,g in protonate( self.Q, molType[0] ):
                if q * aaPerOneCharge < sideChainsNo:
                    yield molType, atomCnt_str, sideChainsNo, q, g

def makeFormulas(fasta, Q, fragType='cz', modifications={}):
    '''Generate all possible fragments given a Roepstorf Scheme [or its generalization].
    '''
    modifications   = standardize(modifications)
    formClass       = {'cz':CZformulator}[fragType](fasta, Q, modifications )
    return formClass

# def genMolecules(fasta, Q, fragmentationScheme='cz', modifications={}, aaPerOneCharge= 5):
#     '''Generate protonated molecules following a given fragmentation scheme.
#     '''
#     IC = isotopeCalculations()
#     precursor, cFrags, zFrags = makeFragments(fasta, fragmentationScheme, modifications)
#     for mol in chain(precursor(),cFrags(),zFrags()):
#         for q,g in protonate( Q, mol['type'] ):
#             if q * aaPerOneCharge < mol['sideChainsNo']:
#                 atomCnt = dict(mol['atomCnt'])
#                 atomCnt['H'] += q + g
#                 monoisotopicMass= IC.getMonoisotopicMass(atomCnt)/float(q)
#                 massMean = IC.getMassMean(atomCnt)/float(q)
#                 massVar  = IC.getMassVar(atomCnt)/float(q**2)
#                 yield ( mol['moleculeType'], q, g, atomCnt, monoisotopicMass, massMean, massVar )

# import pandas as pd
# def pandizeSubstances(precursor, cFrags, zFrags):
#     '''Turns results into a pandas data frame.'''
#     def combineMolecules():
#         for x in chain(precursor(), cFrags(), zFrags()):
#             x['atomCnt']['moleculeType'] = x['moleculeType']
#             yield x['atomCnt']
#     result = pd.DataFrame(combineMolecules()).fillna(0)
#     result[list('CHNOS')] = result[list('CHNOS')].astype(int)
#     idx = ['moleculeType']
#     idx.extend(list('CHNOS'))
#     result = result[idx]
#     return result
