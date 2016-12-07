%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules
from InSilico.spectrumGenerator import insilicoSpectrum, flatten, makeNoise

try:
   import cPickle as pickle
except:
   import pickle

from numpy.random import multinomial
from IsoSpecPy import IsoSpecPy
from Formulator.misc import atomCnt2string

# getting random spectrum:
fasta = substanceP = 'RPKPQQFFGLM'
fasta = ubiquitin = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9

modifications = {}
ionsNo  = 1000000
P       = .999

IS = insilicoSpectrum(fasta, Q, ionsNo, P)



f_charges       = [ float(mol[1]**2) for mol in IS.molecules ]
total_f_charges = sum(f_charges)
probs           = [ q/total_f_charges for q in f_charges ]

moleculeCounts = multinomial( ionsNo, probs )

def genIsotopicEnvelope(isotopologuesNo, atomCnt, jointProb=.999):
    atomCnt_str = atomCnt2string(atomCnt)
    envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
    return envelope

from linearCounter import linearCounter as lcnt
from math import exp

spectrum = lcnt()
for isoCnt, mol in zip(moleculeCounts, IS.molecules):
    if isoCnt>0:
        atomCnt_str = atomCnt2string(mol[3])
        spectrum += lcnt( IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb ) )

isoCnt, mol = moleculeCounts[0], IS.molecules[0]
atomCnt_str = atomCnt2string(mol[3])
envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, .999 )

spectrumDigits = 2

############# it seems that generating these takes more time than other stuff.
%%time
res = lcnt()
for isoCnt, mol in zip(moleculeCounts, IS.molecules):
    atomCnt_str = atomCnt2string(mol[3])
    envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, IS.P )
    for e in envelope.getConfs():
        res[round(e[0], spectrumDigits)] += exp(e[1])
####################################


%%time
spectrum = lcnt()
for isoCnt, mol, prob in zip(moleculeCounts, IS.molecules, probs):
    if isoCnt>0:
        atomCnt_str = atomCnt2string(mol[3])
        envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, IS.P )
        masses, logInts,_ = envelope.getConfsRaw()
        for mz, logint in zip(masses, logInts):
            spectrum[round(mz, spectrumDigits)] += exp(logint)*prob

moleculeCounts = multinomial( ionsNo, spectrum.values() )
moleculeCounts

####################################
import multiprocessing as MP

%%time
def getIsotopes(input):
    isoCnt, mol, prob = input
    res = lcnt()
    if isoCnt>0:
        atomCnt_str = atomCnt2string(mol[3])
        envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, IS.P )
        masses, logInts,_ = envelope.getConfsRaw()
        for mz, logint in zip(masses, logInts):
            res[round(mz, spectrumDigits)] += exp(logint)*prob
    return res

P = MP.Pool(3)
Results = P.map( getIsotopes, zip(moleculeCounts, IS.molecules, probs) )


####################################
%%time
masses = []
logIntensities = []
for isoCnt, mol in zip(moleculeCounts, IS.molecules):
    atomCnt_str = atomCnt2string(mol[3])
    envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, IS.P )
    mass, logInts,_ = envelope.getConfsRaw()
    masses.extend(mass)
    logIntensities.extend(logInts)
###

envelope = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, IS.P )
masses, logIntensities,_ = envelope.getConfsRaw()
# for m, logInt in zip(masses, logIntensities):

lcnt(dict( (round(e[0], spectrumDigits), exp(e[1]))  ))

genIsotopicEnvelope(isoCnt, mol[3], IS.P) for isoCnt, mol in zip(moleculeCounts, IS.molecules)

# %time sample = [ genIsotopicEnvelope(isoCnt, mol[3], IS.P) for isoCnt, mol in zip(moleculeCounts, IS.molecules)]


################################################
############ Checking how big gains when considering no charges.
from Formulator.formulator import makeFragments
from itertools import chain
import numpy as np

jointProb  = .999
prec, c, z = makeFragments(fasta)
prec, c, z = [ list(s()) for s in [prec,c,z] ]

%%time
masses = []
logIntensities = []
for mol in chain(prec,c,z):
    atomCnt_str = atomCnt2string(mol['atomCnt'])
    envelope    = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
    mass, logInts,_ = envelope.getConfsRaw()
    masses.extend(mass)
    logIntensities.extend(logInts)

def getArray(cdata):
    L = len(cdata)
    A = np.empty(L)
    for i, c in enumerate(cdata):
        A[i] = c
    return A

getArray(mass)

for i, el in enumerate(gimme()): my_array[i] = el
np.array(m for m in mass)


for mol in chain(prec,c,z):
    atomCnt_str = atomCnt2string(mol['atomCnt'])
    IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
