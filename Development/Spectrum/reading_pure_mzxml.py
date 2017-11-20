%load_ext autoreload
%autoreload 2
from six.moves import reduce
from MassTodonPy.Spectra.Read import read_mzxml_spectrum,\
                                     read_mzxml_spectrum_faster,\
                                     stack_spectra

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
path = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'


%%timeit
mzxml_spectra = read_mzxml_spectrum_faster(path,
                                           2,
                                           100.0,
                                           True)
spectrum, total_ion_current = reduce(stack_spectra,
                                     mzxml_spectra)

%%timeit
mzxml_spectra = read_mzxml_spectrum(path,
                                    2,
                                    100.0,
                                    True)
spectrum, total_ion_current = reduce(stack_spectra,
                                     mzxml_spectra)
