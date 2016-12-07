import numpy as np
from IsoSpecPy import IsoSpecPy
from math import exp, floor
from collections import Counter
try:
  import cPickle as pickle
except:
  import pickle
from Formulator import makeFragments, protonate
from numpy.random import multinomial
import scipy.stats as ss

def atomCnt2string(atomCnt):
    keys = atomCnt.keys()
    keys.sort()
    return "".join( el+str(atomCnt[el]) for el in keys )

class isotopeCalculator:
    '''A class for isotope calculations.'''
    def __init__(self, massPrecDigits = 3, isoMasses=None, isoProbs=None):
        '''Initiate class with information on isotopes. Calculates basic statistics of isotope frequencies: mean masses and their standard deviations.'''

        if isoMasses==None or isoProbs==None:
            isoMasses, isoProbs = pickle.load(open('data/isotopes.txt', 'rb'))
        self.isoMasses = isoMasses
        self.isoProbs  = isoProbs

        self.elementsMassMean = dict(
            (el, sum( pr*m for pr, m in zip(self.isoProbs[el], self.isoMasses[el]) ) )
            for el in self.isoMasses.keys() )

        self.elementsMassVar  = dict(
            (el, sum( pr*m**2 for pr, m in zip(self.isoProbs[el], self.isoMasses[el])) - self.elementsMassMean[el]**2 )
            for el in self.isoMasses.keys() )

        self.isotopicEnvelopes = {}
        self.massPrecDigits    = massPrecDigits

    def getMonoisotopicMass(self, atomCnt):
        '''Calculate monoisotopic mass of an atom count.'''
        return sum( self.isoMasses[el][0]*elCnt for el, elCnt in atomCnt.items() )

    def getMassMean(self, atomCnt):
        '''Calculate average mass of an atom count.'''
        return sum( self.elementsMassMean[el]*elCnt for el, elCnt in atomCnt.items() )

    def getMassVar(self, atomCnt):
        '''Calculate standard deviation around average mass of an atom count.'''
        return sum( self.elementsMassVar[el]*elCnt for el, elCnt in atomCnt.items() )

    def getOldEnvelope(self, atomCnt_str):
        masses, probs = self.isotopicEnvelopes[atomCnt_str]
        return masses.copy(), probs.copy()

    def isoEnvelope(self, atomCnt, jointProb, q=0, g=0):
        '''Get an isotopic envelope consisting of a numpy array of masses and numpy array of probabilities.'''

        atomCnt_str = atomCnt2string(atomCnt)
        if atomCnt_str in self.isotopicEnvelopes:
            masses, probs = self.getOldEnvelope(atomCnt_str)
        else:
            counts = []
            isotope_masses = []
            isotope_probs  = []
            for el, cnt in atomCnt.items():
                counts.append(cnt)
                isotope_masses.append(self.isoMasses[el])
                isotope_probs.append(self.isoProbs[el])
            envelope = IsoSpecPy.IsoSpec( counts, isotope_masses, isotope_probs, jointProb )

            masses, logprobs, _ = envelope.getConfsRaw()
            aggregator = Counter()
            for mass, logprob in zip(masses, logprobs):
                aggregator[ round(mass, self.massPrecDigits) ] += exp(logprob)

            masses = np.array(aggregator.keys())
            probs  = np.empty(len(masses))
            for i, mass in enumerate(masses):
                probs[i] = aggregator[mass]

            self.isotopicEnvelopes[atomCnt_str] = ( masses.copy(), probs.copy() ) # memoization

        masses = np.around( (masses + g + q)/q, decimals=self.massPrecDigits )
        return masses, probs
    #TODO add a version that performs all possible calculations.
    #TODO add a version that uses precalculated spectra for some substances like proteins/metabolites/so on .. so on.. This would save massively time for generation.

    def randomSpectrum(self, fasta, Q, ionsNo, fragScheme='cz', aaPerOneCharge=5, jointProb=.999, scale =.01, modifications={} ):
        '''Get random spectrum following a heuristical data generation process.'''

        averageSpectrum     = Counter()
        chargesSquaredSum   = 0.0
        precs, cfrags, zfrags = makeFragments( fasta, fragScheme, modifications )
        for mol in chain( precs(), cfrags(), zfrags() ):
            for q,g in protonate( Q, mol['type'] ):
                if q * aaPerOneCharge < mol['sideChainsNo']:
                    chargesSquaredSum += q**2
                    masses, probs = self.isoEnvelope( mol['atomCnt'], jointProb, q, g )
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

        spectrum = np.around(spectrum, self.massPrecDigits)
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
