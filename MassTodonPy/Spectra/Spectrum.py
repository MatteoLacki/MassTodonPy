"""Deal with experimental spectrum."""

from MassTodonPy.Data.Constants import eps
from MassTodonPy.Spectra.Read import read_spectrum
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum

class Spectrum(object):
    """Prepare experimental spectrum.

    Parameters
    ----------
    spectrum : string or tuple of numpy arrays or ExperimentalSpectrum
        The path can end up with extension:
            *.txt, for spectra saved in tab separated format.
            *.mzXml, for spectra saved with mxXml format.
        The tuple and experimental spectrum consist of two numpy arrays,
        containing m over z ratios and respective intensities.
    mz_precision : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    minimal_intensity : float
        Experimental peaks with lower height will be trimmed.
    percent_top_peaks : float
        Percentage of the heighest peaks in the spectrum to be included.

    """

    def __init__(self,
                 spectrum='',
                 mz_precision=2,
                 minimal_intensity=eps,
                 percent_top_peaks=1.0):
        """Initialize the Spectrum."""
        self.mz_precision = mz_precision
        if isinstance(spectrum, str):  # read in spectrum from file
            self.spectrum = read_spectrum(spectrum, mz_precision, eps)
        else:  # use the provided spectrum
            if isinstance(spectrum, tuple):
                self.spectrum = ExperimentalSpectrum(*spectrum)
            elif isinstance(spectrum, ExperimentalSpectrum):
                self.spectrum = spectrum
            self.spectrum.trim_intensity(eps)  # first: cut zero intensity m/z
            self.spectrum.round_masses(mz_precision)

        self.total_intensity = self.spectrum.total_intensity()
        self.left_over_spectrum = 0  # outside the analysis
        if minimal_intensity > eps:
            self.left_over_spectrum = self.spectrum.split_measure(minimal_intensity)
        if 0 < percent_top_peaks < 1:
            cut_off = self.spectrum.get_P_set_cut_off(percent_top_peaks)
            self.left_over_spectrum += self.spectrum.split_measure(cut_off)
        self.intensity_after_trim = self.spectrum.total_intensity()

    def __iter__(self):
        """Iterate over the spectrum."""
        return self.spectrum.__iter__()

    def _low_spectrum(self):
        """Iterate over left-over spectrum."""
        return self.left_over_spectrum.__iter__()

    def __getitem__(self, mz_interval):
        """Filter part of the spectrum."""
        return self.spectrum.__getitem__(mz_interval)
