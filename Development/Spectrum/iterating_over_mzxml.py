from pyteomics import mzxml
from collections import defaultdict, Counter

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
path = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'

entries = defaultdict(Counter)

with mzxml.read(path) as reader:
    for spectrum in reader:
        mzs = spectrum['m/z array']
        intensities = spectrum['intensity array']
        entries['msLevel'][spectrum['msLevel']] += 1
        precursorMz = spectrum['precursorMz'][0]
        entries['precursorMz'][precursorMz['precursorMz']] += 1
        entries['windowWideness'][precursorMz['windowWideness']] += 1
        entries['centroided'][spectrum['centroided']] += 1
spectrum
