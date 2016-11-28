from IsoSpecPy import IsoSpecPy
from Formulator.isotopeCalculator import IsotopeCalculations
from numpy.random import multinomial
from math import exp
from Formulator.formulator import genMolecules
from Formulator.misc import atomCnt2string
from itertools import chain
import scipy.stats as ss

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
        '''Simulates a mass spectrumin silico.'''
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
        getMols = ( genIsotopicEnvelope(isoCnt, mol[3], self.P)
            for isoCnt, mol in zip(moleculeCounts, self.molecules) )
        result = []
        for m,n in chain(*getMols):
            result.extend(ss.norm.rvs(loc=m, size=n) )
        return result
