%load_ext autoreload
%autoreload 2

from MassTodonPy.Spectra.Measure import Measure
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Spectra.Read import _read_mz

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
path = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'

%%time
spectra = _read_mz(path, 3, 'mzxml', 100.0, sum_intensity=True)

spectrum = next(spectra)
spectrum
x = "aBc"
x.lower()
x


x = 'a'

if x is 'c' or 'a':
    print('Hurra')
else:
    print('Buuu')
'Buuu'[1:]
