from IsoSpecPy import IsoSpecPy
from Formulator.isotopeCalculator import IsotopeCalculations
from numpy.random import multinomial
from math import exp, floor
from Formulator.formulator import genMolecules
from Formulator.misc import atomCnt2string
from itertools import chain
import scipy.stats as ss
from operator import itemgetter


def genIsotopicEnvelope(isotopologuesNo, atomCnt, jointProb=.999):
    atomCnt_str = atomCnt2string(atomCnt)
    envelope = IsoSpecPy.IsoSpec.IsoFromFormula( atomCnt_str, jointProb )
    masses = []
    probs = []
    for x in envelope.getConfs():
        masses.append(x[0])
        probs.append(exp(x[1]))
    counts = multinomial( isotopologuesNo, probs )
    return zip(masses, counts)

class insilicoSpectrum:
    def __init__(self, fasta, Q, ionsNo, jointProb=.999, fragScheme='cz', isoMasses=None, isoProbs=None, modifications={}):
        '''Simulate a mass spectrum in silico or in blanco.'''
        self.P = jointProb
        self.fasta = fasta
        self.Q = Q
        self.molecules = list( genMolecules(fasta, Q, fragScheme, modifications) )

    def rvs(self, ionsNo):
        '''Generate random mass spectrum.'''
        f_charges = [ float(mol[1]**2) for mol in self.molecules ] # Probability of being chosen proportional to the square of number of charges.
        total_f_charges = sum(f_charges)
        probs = [ q/total_f_charges for q in f_charges ]
        moleculeCounts = multinomial( ionsNo, probs )
        return [
            list(genIsotopicEnvelope(isoCnt, mol[3], self.P))
            for isoCnt, mol in zip(moleculeCounts, self.molecules)
        ], probs

def flatten( massSpectra ):
    result = []
    for sp in massSpectra:
        result.extend(sp)
    return result

def makeNoise( MassSpectrum, percentPeaks = .2 ):
    '''Produce noise peaks using a strategy that is totally atheoretic.'''
    M = max(MassSpectrum, key=itemgetter(0))[0]
    Imean   = sum(i for m,i in MassSpectrum)/len(MassSpectrum)
    size    = int(floor(len(MassSpectrum)*percentPeaks))
    masses = ss.uniform.rvs(    loc     = .0,
                                scale   = 1.1*M,
                                size    = size  )
    intensities = ss.poisson.rvs(mu=Imean, size=size )
    return list(zip(masses, intensities))
