"""Deal with experimental spectrum."""
import csv
import numpy as np

from MassTodonPy.Data.Constants import eps, infinity
from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Spectra.Read import read_spectrum
from MassTodonPy.Measure.Measure import Measure


class Spectrum(Measure):
    """Prepare experimental spectrum.

    Parameters
    ----------
    spectrum : string, tuple, list, or Spectrum
        The path can end up with extension:
            *.txt, for spectra saved in tab separated format.
            *.mzXml, for spectra saved with mxXml format.
        The tuple or list consist of two numpy arrays,
        one with m over z ratios and the other with intensities.
    mz_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    min_intensity : float
        Experimental peaks with lower height will be trimmed.
    percent_top_peaks : float
        Percentage of the heighest peaks in the spectrum to be included.

    """
    def __init__(self,
                 mz=np.array([]), 
                 intensity=np.array([]),
                 spectrum='',
                 mz_digits=2,
                 min_intensity=eps,
                 percent_top_peaks=1.0):
        """Initialize the Spectrum."""
        self._store_names = ('mz', 'intensity')
        self.mz_digits = mz_digits
        if isinstance(spectrum, str) and spectrum: 
            spectrum = read_spectrum(spectrum, mz_digits, eps)
            self.mz = spectrum.atoms
            self.intensity = spectrum.masses
        elif isinstance(spectrum, (Measure, Spectrum)):  # copy 'constructor'
            self.mz = spectrum.atoms.copy()
            self.intensity = spectrum.masses.copy()
        else:
            if isinstance(spectrum, (tuple, list)):
                mz, intensity = spectrum
            elif not spectrum:  # potentially an empty spectrum
                mz = mz
                intensity = intensity            
            mz = np.array(mz)
            intensity = np.array(intensity)
            assert all(intensity >= 0), "intensity must be non-negative"
            assert len(mz) == len(intensity)
            self.mz = mz
            self.intensity = intensity
            self.trim_intensity(eps)  # intensity > 0
            self.round_mz(mz_digits)
        self.low_spectrum = 0
        if min_intensity > eps:
            self.low_spectrum += self.split_measure(min_intensity)
        if 0 < percent_top_peaks < 1:
            cut_off = self.get_P_set_cut_off(percent_top_peaks)
            self.low_spectrum += self.split_measure(cut_off)

    def split_measure(self, cut_off):
        """Split measure into two according to the cut off on masses.

        Retain the measure with masses greater or equal to the cut off.
        Parameters
        ----------
        cut_off : float
        Returns
        ----------
        other : Measure
            A measure with masses strictly below the cut off.
        """
        mz = self.mz[self.intensity < cut_off]
        intensity = self.intensity[self.intensity < cut_off]
        return Spectrum(mz, intensity, '', self.mz_digits, eps, 1.0)

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def intensity(self):
        """Get intensities."""
        return self.masses

    @intensity.setter
    def intensity(self, intensity):
        """Set intensities."""
        self.masses = intensity

    def round_mz(self, precision=infinity):
        """Round the atoms of the measure to a given precision.

        Parameters
        ----------
        precision : integer
            The number of digits after which the atoms' masses get rounded.
            E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
            Defaults to 'inf', which prevents any rounding.

        """
        self.round_atoms(precision)

    def total_intensity(self):
        return self.intensity.sum()

    def trim_intensity(self, cut_off):
        """Trim intensities below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        self.trim(cut_off)
