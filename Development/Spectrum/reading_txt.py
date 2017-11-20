%load_ext autoreload
%autoreload 2
from six.moves import reduce
from MassTodonPy.Spectra.Read import read_txt_spectrum,\
                                     stack_spectra

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_txt = spectrum_path + 'subP_spectrum.txt'


read_txt_spectrum(spectrum_txt, 2, 100.0, True)
