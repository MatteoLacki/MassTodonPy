from IsoSpecPy import IsoSpecPy
from Formulator.isotopeCalculator import IsotopeCalculations
from numpy.random import multinomial

class insilicoSpectrum:
    def __init__(self, fasta, ionsNo, fragScheme='cz', jointProb=.999, isoMasses=None, isoProbs=None):
        '''Simulates in silico mass spectrum.'''
        self.N = ionsNo
        self.P = jointProb
        self.fasta = fasta

    def rvs():
        pass
