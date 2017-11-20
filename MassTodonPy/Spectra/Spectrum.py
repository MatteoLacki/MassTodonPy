from six.moves import reduce
from MassTodonPy.Spectra.Operations import parse_path,\
                                           round_n_trim,\
                                           stack_spectra,\
                                           get_minimal_intensity
from MassTodonPy.Spectra.Read import read_mzxml_spectrum,\
                                     read_mzxml_spectrum_faster,\
                                     read_txt_spectrum
from MassTodonPy.Spectra.Spectra import ExperimentalSpectrum as ExpSpec
from MassTodonPy.Data.Constants import eps


class Spectrum(object):
    """
    Store experimental spectrum.
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
        assert spectrum is not "", "Provide a spectrum to analyze!"
        assert isinstance(spectrum, (tuple, str, ExpSpec)),\
            "spectrum is not a string,\
            nor a tuple, not even an ExperimentalSpectrum."
        self.spectrum = spectrum
        self.mz_precision_digits = mz_precision_digits
        self.minimal_intensity = minimal_intensity
        self.percent_top_peaks = percent_top_peaks
        self._faster = _faster
        if isinstance(spectrum, str):
            self.read_and_summarize_spectrum(eps)
        else:
            self.round_n_trim_spectrum(eps)
        if 0 < percent_top_peaks < 1:
            self.minimal_intensity = get_minimal_intensity(
                *self.spectrum,
                retained_percentage=percent_top_peaks)
        self.intensity_after_trim = sum(intensity for mz, intensity in
                                        self.high_spectrum())

    def read_and_summarize_spectrum(self, minimal_intensity):
        file_path, file_name, file_ext = parse_path(self.spectrum)
        file_ext = file_ext.lower()
        if file_ext is '':
            file_ext = '.txt'
        assert file_ext in ('.txt', '.mzxml'), "Unsupported extension."
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
            raise NotImplementedError
        else:
            raise NotImplementedError

    def round_n_trim_spectrum(self, minimal_intensity):
        self.total_intensity = sum(self.spectrum[1])
        self.spectrum = ExpSpec(*round_n_trim(
            *self.spectrum,
            support_precision=self.mz_precision_digits,
            value_cut_off=minimal_intensity))

    def high_spectrum(self):
        for mz, intensity in zip(*self.spectrum):
            if intensity >= self.minimal_intensity:
                yield mz, intensity

    def low_spectrum(self):
        for mz, intensity in zip(*self.spectrum):
            if intensity < self.minimal_intensity:
                yield mz, intensity
