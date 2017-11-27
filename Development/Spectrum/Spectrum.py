%load_ext autoreload
%autoreload 2

from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.Data.Constants import eps, infinity

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_mzxml = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'

%%timeit
ms = Spectrum(spectrum=spectrum_mzxml, mz_precision=2,
              minimal_intensity=100.0, percent_top_peaks=.95)

ms.spectrum

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_txt = spectrum_path + 'subP_spectrum.txt'

%%time
ms = Spectrum(spectrum=spectrum_txt, mz_precision=3,
              minimal_intensity=100.0, percent_top_peaks=.95)
ms.spectrum

%%time
ms = Spectrum(spectrum=spectrum_txt, mz_precision=2,
              minimal_intensity=100.0, percent_top_peaks=1.0)
ms.spectrum

from MassTodonPy.Data.get_dataset import get_dataset
mol = get_dataset('substanceP')

ms = Spectrum(spectrum=mol.spectrum, mz_precision=2,
              minimal_intensity=100.0, percent_top_peaks=.95)
ms.spectrum

# mol = get_dataset('ubiquitin')
#
# ms = Spectrum(spectrum=mol.spectrum, mz_precision=2,
#               minimal_intensity=100.0, percent_top_peaks=.95)
#
# list(ms.high_spectrum())


if None:
    print('Ha')
else:
    print('Be')
