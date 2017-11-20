from six.moves import reduce
from MassTodonPy.Spectra.operations import parse_path, round_n_trim
from MassTodonPy.Spectra.Read import stack_spectra,\
                                     read_mzxml_spectrum,\
                                     read_mzxml_spectrum_faster,\
                                     read_txt_spectrum
from MassTodonPy.Spectra.Spectra import ExperimentalSpectrum as ExpSpec


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
    precision_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    """

    def __init__(self,
                 spectrum='',
                 intensity_cut_off=0.0,
                 precision_digits=2,
                 _faster=False):
        assert spectrum is not "", "Provide a spectrum to analyze!"
        assert isinstance(spectrum, (tuple, str, ExpSpec)),\
            "spectrum is not a string,\
            nor a tuple, not even an ExperimentalSpectrum."
        self.intensity_cut_off = intensity_cut_off
        self.precision_digits = precision_digits
        self.spectrum = spectrum
        self._faster = _faster
        if isinstance(spectrum, str):
            self.read_and_summarize_spectrum()
        else:
            self.round_n_trim()

    def read_and_summarize_spectrum(self):
        file_path, file_name, file_ext = parse_path(self.spectrum)
        file_ext = file_ext.lower()
        if file_ext is '':
            file_ext = '.txt'
        assert file_ext in ('.txt', '.mzxml'), "Unsupported extension."
        if file_ext == '.txt':
            self.spectrum, self.total_ion_current = read_txt_spectrum(
                self.spectrum,
                self.precision_digits,
                self.intensity_cut_off,
                sum_intensity=True)
        elif file_ext == '.mzxml':
            if self._faster:
                mzxml_spectra = read_mzxml_spectrum_faster
            else:
                mzxml_spectra = read_mzxml_spectrum
            mzxml_spectra = mzxml_spectra(self.spectrum,
                                          self.precision_digits,
                                          self.intensity_cut_off,
                                          sum_intensity=True)
            self.spectrum, self.total_ion_current = reduce(stack_spectra,
                                                           mzxml_spectra)
        elif file_ext == '.mzml':
            raise NotImplementedError
        else:
            raise NotImplementedError

    def round_n_trim(self):
        self.total_ion_current = sum(self.spectrum[1])
        self.spectrum = ExpSpec(*round_n_trim(*self.spectrum,
                                              self.precision_digits,
                                              self.intensity_cut_off))
