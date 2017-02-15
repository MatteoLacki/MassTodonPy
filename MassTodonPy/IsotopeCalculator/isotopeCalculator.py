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

from IsoSpecPy      import IsoSpecPy
from formulaParser  import formulaParser
from math           import exp, floor, fsum
from collections    import Counter, defaultdict
try:
  import cPickle as pickle
except:
  import pickle
from numpy.random   import multinomial
import scipy.stats  as ss
import numpy        as np


def cdata2numpyarray(x):
    '''Turn c-data into a numpy array.'''
    res = np.empty(len(x))
    for i in xrange(len(x)):
        res[i] = x[i]
    return res

def agg_spec_proper(masses, probs, digits=2):
    '''Aggregate values with the same keys.'''
    lists = defaultdict(list)
    for mass, prob in zip(masses.round(digits), probs):
        lists[mass].append(prob)
    newMasses = np.array(lists.keys())
    newProbs  = np.empty(len(newMasses))
    for prob, mass in zip(np.nditer(newProbs,op_flags=['readwrite']), newMasses):
        prob[...] = fsum(lists[mass])
    return newMasses, newProbs

def aggregate( keys, values ):
    '''Aggregate values with the same keys.'''
    uniqueKeys, indices = np.unique( keys, return_inverse=True)
    return uniqueKeys, np.bincount( indices, weights=values )

class isotopeCalculator:
    '''A class for isotope calculations.'''
    def __init__(   self,
                    jP = .999,
                    precDigits= 2,
                    isoMasses = None,
                    isoProbs  = None  ):
        '''Initiate class with information on isotopes. Calculates basic statistics of isotope frequencies: mean masses and their standard deviations.'''

        if isoMasses==None or isoProbs==None:
            isoMasses, isoProbs = pickle.load(open('data/isotopes.txt', 'rb'))
        self.isoMasses = isoMasses
        self.isoProbs  = isoProbs
        self.jP = jP

        self.elementsMassMean = dict(
            (el, sum( pr*m for pr, m in zip(self.isoProbs[el], self.isoMasses[el]) ) )
            for el in self.isoMasses.keys() )

        self.elementsMassVar  = dict(
            (el, sum( pr*m**2 for pr, m in zip(self.isoProbs[el], self.isoMasses[el])) - self.elementsMassMean[el]**2 )
            for el in self.isoMasses.keys() )

        self.isotopicEnvelopes = {}
        self.precDigits = precDigits
        self.formParser = formulaParser()

    def getMonoisotopicMass(self, atomCnt):
        '''Calculate monoisotopic mass of an atom count.'''
        return sum( self.isoMasses[el][0]*elCnt for el, elCnt in atomCnt.items() )

    def getMassMean(self, atomCnt):
        '''Calculate average mass of an atom count.'''
        return sum( self.elementsMassMean[el]*elCnt for el, elCnt in atomCnt.items() )

    def getMassVar(self, atomCnt):
        '''Calculate mass variance of an atom count.'''
        return sum( self.elementsMassVar[el]*elCnt for el, elCnt in atomCnt.items() )

    def getSummary(self, atomCnt_str):
        atomCnt = self.formParser.parse(atomCnt_str)
        return (    self.getMonoisotopicMass(atomCnt),
                    self.getMassMean(atomCnt),
                    self.getMassVar(atomCnt)    )

    def getOldEnvelope(self, atomCnt_str, jP, precDigits):
        masses, probs = self.isotopicEnvelopes[(atomCnt_str, jP, precDigits)]
        return masses.copy(), probs.copy()

    def getNewEnvelope(self, atomCnt_str, jP, precDigits):
        counts          = []
        isotope_masses  = []
        isotope_probs   = []
        atomCnt = self.formParser.parse(atomCnt_str)
        for el, cnt in atomCnt.items():
            counts.append(cnt)
            isotope_masses.append(self.isoMasses[el])
            isotope_probs.append(self.isoProbs[el])
        envelope = IsoSpecPy.IsoSpec( counts, isotope_masses, isotope_probs, jP )
        masses, logprobs, _ = envelope.getConfsRaw()
        masses  = cdata2numpyarray(masses)
        probs   = np.exp(cdata2numpyarray(logprobs))
        masses, probs = agg_spec_proper(masses, probs, precDigits)
        # memoization
        self.isotopicEnvelopes[ (atomCnt_str, jP, precDigits) ] = ( masses, probs )
        return masses.copy(), probs.copy()

    def getEnvelope(self, atomCnt_str, jP, precDigits):
        if (atomCnt_str, jP, precDigits) in self.isotopicEnvelopes:
            masses, probs = self.getOldEnvelope(atomCnt_str,jP,precDigits)
        else:
            masses, probs = self.getNewEnvelope(atomCnt_str,jP,precDigits)
        return masses, probs

    def isoEnvelope(self, atomCnt_str, jP=None, q=0, g=0, precDigits=None):
        '''Get an isotopic envelope consisting of a numpy array of masses and numpy array of probabilities.'''

        if jP is None:
            jP = self.jP
        if precDigits is None:
            precDigits = self.precDigits
        masses, probs  = self.getEnvelope(atomCnt_str, jP, precDigits)
        if q is not 0:
            masses = np.around( (masses + g + q)/q, decimals=precDigits )
        masses, probs = aggregate(masses, probs)
        return masses, probs


    def randomFragmentationExperiment(self, fasta, Q, ionsNo, fragmentator, aaPerOneCharge=5, jointProb=.999, scale =.01, modifications={} ):
        '''Get random spectrum of a fragmentation experiment.'''

        averageSpectrum     = Counter()
        chargesSquaredSum   = 0.0

        for molType, atomCnt_str, sideChainsNo, q, g in fragmentator.makeMolecules(aaPerOneCharge):
            chargesSquaredSum += q**2
            masses, probs = self.isoEnvelope( atomCnt_str, jointProb, q, g )
            for mass, prob in zip(masses,probs):
                averageSpectrum[mass] += prob * q**2

        masses = averageSpectrum.keys()
        probs  = np.empty(len(masses))
        for i, m in enumerate(masses):
            probs[i] = averageSpectrum[m]
        probs = probs/chargesSquaredSum
        ionsPerMass = multinomial(ionsNo, probs)
        spectrum    = np.empty(ionsNo)
        i = 0
        for mass, ionCnt in zip(masses, ionsPerMass):
            if ionCnt > 0:
                for mz in np.random.normal( loc=mass, scale=.01, size=ionCnt ):
                    spectrum[i] = mz
                    i += 1

        spectrum = np.around(spectrum, self.precDigits)
        spectrum = Counter(spectrum)
        masses   = np.array(spectrum.keys())
        intensities = np.empty(len(masses))
        for i, m in enumerate(masses):
            intensities[i] = spectrum[m]
        return masses, intensities

    def randomPrecursors(self, fasta, Q, ionsNo, fragmentator, aaPerOneCharge=5, jointProb=.999, scale =.01, modifications={} ):
        '''Get random spectrum of a fragmentation experiment.'''

        averageSpectrum     = Counter()
        chargesSquaredSum   = 0.0

        for molType, atomCnt_str, sideChainsNo, q, g in fragmentator.makeMolecules(aaPerOneCharge):
            if molType == 'precursor':
                chargesSquaredSum += q**2
                masses, probs = self.isoEnvelope( atomCnt_str, jointProb, q, g )
                for mass, prob in zip(masses,probs):
                    averageSpectrum[mass] += prob * q**2

        masses = averageSpectrum.keys()
        probs  = np.empty(len(masses))
        for i, m in enumerate(masses):
            probs[i] = averageSpectrum[m]
        probs = probs/chargesSquaredSum
        ionsPerMass = multinomial(ionsNo, probs)
        spectrum    = np.empty(ionsNo)
        i = 0
        for mass, ionCnt in zip(masses, ionsPerMass):
            if ionCnt > 0:
                for mz in np.random.normal( loc=mass, scale=.01, size=ionCnt ):
                    spectrum[i] = mz
                    i += 1

        spectrum = np.around(spectrum,self.precDigits)
        spectrum = Counter(spectrum)
        masses   = np.array(spectrum.keys())
        intensities = np.empty(len(masses))
        for i, m in enumerate(masses):
            intensities[i] = spectrum[m]
        return masses, intensities

    def addNoise(self, masses, intensities, percentPeaks = .2):
        '''Produce noise peaks using a strategy that is totally atheoretic.'''
        M       = max(masses)
        Imean   = intensities.mean()
        size    = int(floor(len(masses)*percentPeaks))
        noise_masses = ss.uniform.rvs(
            loc     = .0,
            scale   = 1.1*M,
            size    = size  )
        noise_intensities = ss.poisson.rvs(mu=Imean, size=size )
        return noise_masses, noise_intensities

    def randomSpectrum(self, atomCnt_str, ionsNo, digits=2, jointProb=.9999, Q=0, sigma=None):
        '''Generate a random spectrum.'''
        #TODO this should be all in all replaced by Michał's software that uses online data generation.
        # atomCnt_str     = atomCnt2string(atomCnt)
        masses, probs   = self.getEnvelope(atomCnt, jointProb, digits)
        counts  = multinomial(ionsNo, probs)
        masses  = masses[ counts > 0 ]
        counts  = counts[ counts > 0 ]
        if Q > 0:
            masses  = (masses + Q)/Q
        noise_masses = np.empty(ionsNo)
        if not sigma==None:
            sigma   = .01
            cnt     = 0
            for mass, c in zip(masses, counts):
                noise_masses[cnt:(cnt+c)] = np.random.normal(mass, sigma, c)
                cnt += c
            noise_masses = noise_masses.round(digits)
            masses, counts = np.unique(noise_masses, return_counts=True)
        return masses, counts
