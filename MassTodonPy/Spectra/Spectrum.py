"""Class to deal with experimental spectrum."""

from six.moves import reduce

from MassTodonPy.Data.Constants import eps
from MassTodonPy.Spectra.Operations import get_minimal_intensity
from MassTodonPy.Spectra.Operations import parse_path
from MassTodonPy.Spectra.Operations import round_n_trim
from MassTodonPy.Spectra.Operations import stack_spectra
from MassTodonPy.Spectra.Read import read_mzml_spectrum
from MassTodonPy.Spectra.Read import read_mzxml_spectrum
from MassTodonPy.Spectra.Read import read_mzxml_spectrum_faster
from MassTodonPy.Spectra.Read import read_txt_spectrum
from MassTodonPy.Spectra.Spectra import ExperimentalSpectrum as ExpSpec


def successor(X, value):
    u"""Find minimal i s.t. X[i] ≥ value."""
    L = 0
    R = len(X)-1
    while L + 1 < R:
        I = (L+R)//2
        if X[I] is value:
            return I
        elif X[I] > value:
            R = I
        else:
            L = I
    if X[I] > value:
        return L
    else:
        return R


class Spectrum(object):
    """Store experimental spectrum.

    Parameters
    ----------
    spectrum : path string or tuple of numpy arrays or ExperimentalSpectrum
        The path can end up with extension:
            *.txt, for spectra saved in tab separated format.
            *.mzXml, for spectra saved with mxXml format.
        The tuple and experimental spectrum consist of two numpy arrays,
        containing m over z ratios and respective intensities.
    mz_precision_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    minimal_intensity : float
        Experimental peaks with lower height will be trimmed.
    percent_top_peaks : float
        Percentage of the heighest peaks in the spectrum to be included.
    _faster : bool
        Use the direct Xml tree parser.

    """

    def __init__(self,
                 spectrum='',
                 mz_precision_digits=2,
                 minimal_intensity=eps,
                 percent_top_peaks=1.0,
                 _faster=False):
        """Initialize the Spectrum."""
        self.spectrum = spectrum
        self.mz_precision_digits = mz_precision_digits
        self.minimal_intensity = minimal_intensity
        self.percent_top_peaks = percent_top_peaks
        self._faster = _faster
        if isinstance(spectrum, str):
            self.read_and_summarize(eps)
        else:
            self.round_n_trim(eps)
        if 0 < percent_top_peaks < 1:
            self.minimal_intensity = get_minimal_intensity(
                *self.spectrum,
                retained_percentage=percent_top_peaks)
        self.get_high_spectrum()
        self.intensity_after_trim = sum(self.high_spectrum.intensity)

    def read_and_summarize(self,
                           minimal_intensity=eps):
        """Read and summarize the spectrum.

        Read from either .txt or .mzXML format.
        """
        assert self.spectrum is not '', "Provide a spectrum to analyze!"
        file_path, file_name, file_ext = parse_path(self.spectrum)
        file_ext = file_ext.lower()
        if file_ext is '':
            file_ext = '.txt'
        if file_ext == '.txt':
            self.spectrum, self.total_intensity = read_txt_spectrum(
                self.spectrum,
                self.mz_precision_digits,
                minimal_intensity,
                sum_intensity=True)
        elif file_ext == '.mzxml':
            if self._faster:
                mzxml_spectra = read_mzxml_spectrum_faster
            else:
                mzxml_spectra = read_mzxml_spectrum
            mzxml_spectra = mzxml_spectra(self.spectrum,
                                          self.mz_precision_digits,
                                          minimal_intensity,
                                          sum_intensity=True)
            self.spectrum, self.total_intensity = reduce(stack_spectra,
                                                         mzxml_spectra)
        elif file_ext == '.mzml':
            mzml_spectra = read_mzml_spectrum(self.spectrum,
                                              self.mz_precision_digits,
                                              minimal_intensity,
                                              sum_intensity=True)
            self.spectrum, self.total_intensity = reduce(stack_spectra,
                                                         mzml_spectra)
        else:
            raise NotImplementedError

    def round_n_trim(self, minimal_intensity):
        """Round and trim spectrum.

        Rount the spectrum to mz_precision_digits and cut peaks below the
        minimal intensity.
        """
        self.total_intensity = sum(self.spectrum[1])
        self.spectrum = ExpSpec(*round_n_trim(
            *self.spectrum,
            support_precision=self.mz_precision_digits,
            value_cut_off=minimal_intensity))

    def get_high_spectrum(self):
        """Get spectrum higher than minimal intensity."""
        mzs, intensities = self.spectrum
        self.high_spectrum = ExpSpec(
            mz=mzs[intensities >= self.minimal_intensity],
            intensity=intensities[intensities >= self.minimal_intensity])

    def _high_spectrum(self):
        """Iterate over spectrum above minimal intensity."""
        for mz, intensity in zip(*self.spectrum):
            if intensity >= self.minimal_intensity:
                yield mz, intensity

    def __iter__(self):
        """Iterate over spectrum above minimal intensity."""
        return self._high_spectrum()

    def _low_spectrum(self):
        """Iterate over spectrum below minimal intensity."""
        for mz, intensity in zip(*self.spectrum):
            if intensity < self.minimal_intensity:
                yield mz, intensity

    def __getitem__(self, mz_interval):
        """Filter part of the spectrum."""
        try:
            mz_min, mz_max = mz_interval
            idx = successor(self.spectrum, mz_min)
            while idx <= mz_max:
                pass
                # yield self.
        except TypeError:
            print("Key must be a tuple with minimal and maximal m/z.")
