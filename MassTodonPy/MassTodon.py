from Formulator         import makeFragments
from IsotopeCalculator  import isotopeCalculator
from PeakPicker         import peakPicker

class MassTodon():
    def __init__( self, fasta, precursorCharge, modifications={}, massPrecDigits = 3, isoMasses=None, isoProbs=None ):
        self.fasta  = fasta
        self.Q      = precursorCharge
        self.isoCalc= isotopeCalculator(massPrecDigits, isoMasses, isoProbs)

    def randomSpectrum(
            self,
            ionsNo,
            fragScheme      = 'cz',
            aaPerOneCharge  = 5,
            jointProb       = .999,
            scale           = .01,
            percentPeaks    =.2     ):

        masses, intensities = self.isoCalc.randomSpectrum(
            self.fasta,
            self.Q,
            ionsNo,
            fragScheme='cz',
            aaPerOneCharge=5,
            jointProb=.999,
            scale =.01          )

        noise_masses, noise_intensities = self.isoCalc.addNoise(
            masses,
            intensities,
            percentPeaks        )

        return masses, intensities, noise_masses, noise_intensities

    def setMassSpectrum(self, massSpectrum):
        self.massSpectrum = massSpectrum

    def pickPeaks(self):
        self.peakPicker = PeakPicker()
