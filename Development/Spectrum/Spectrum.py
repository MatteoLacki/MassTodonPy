%load_ext autoreload
%autoreload 2


from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.Data.Constants import eps

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_mzxml = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'

%%time
ms = Spectrum(spectrum=spectrum_mzxml, mz_precision_digits=2,
              minimal_intensity=100.0, percent_top_peaks=.95)
ms.minimal_intensity
ms.spectrum
ms.high_spectrum

list(ms.high_spectrum())

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_txt = spectrum_path + 'subP_spectrum.txt'

%%time
ms = Spectrum(spectrum=spectrum_txt, mz_precision_digits=2,
              minimal_intensity=100.0, percent_top_peaks=.95)
ms.minimal_intensity

%%time
ms = Spectrum(spectrum=spectrum_txt, mz_precision_digits=2,
              minimal_intensity=100.0, percent_top_peaks=1.0)
ms.minimal_intensity

from MassTodonPy.Data.get_dataset import get_dataset
mol = get_dataset('substanceP')

ms = Spectrum(spectrum=mol.spectrum, mz_precision_digits=2,
              minimal_intensity=100.0, percent_top_peaks=.95)

list(ms.high_spectrum())



mol = get_dataset('ubiquitin')

ms = Spectrum(spectrum=mol.spectrum, mz_precision_digits=2,
              minimal_intensity=100.0, percent_top_peaks=.95)

list(ms.high_spectrum())
