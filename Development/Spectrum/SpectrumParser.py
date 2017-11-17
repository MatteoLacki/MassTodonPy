%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Spectra.SpectrumParser import \
    read_n_preprocess_spectrum

subP = get_dataset('substanceP')
spectrum_from_spectrum = read_n_preprocess_spectrum(subP.spectrum,
                                                    100., 1.0, 2)

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'

spectrum_txt = spectrum_path + 'subP_spectrum.txt'
spectrum_from_txt = read_n_preprocess_spectrum(spectrum_txt,
                                               100., 1.0, 2)

spectrum_mzxml = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'
spectrum_from_mzxml = read_n_preprocess_spectrum(spectrum_mzxml,
                                                 100., 1.0, 2)
