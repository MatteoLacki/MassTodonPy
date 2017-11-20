%load_ext autoreload
%autoreload 2

from MassTodonPy.Spectra.Spectrum import Spectrum

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
path = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'

%%timeit
ms = Spectrum(path, 100.0, 2)
ms.spectrum
