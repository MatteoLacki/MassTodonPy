from pyteomics          import mzxml
from IsotopeCalculator  import aggregate
import numpy as np
import os

def merge_runs(spec1, spec2):
    mz = np.concatenate((spec1[0], spec2[0]))
    I  = np.concatenate((spec1[1], spec2[1]))
    return aggregate(mz, I)

def get_spectra(path, cutOff=100, digits=2):
    with mzxml.read(path) as reader:
        for spectrum in reader:
            mz          = spectrum['m/z array']
            intensity   = spectrum['intensity array']
            mz          = mz[intensity > cutOff]
            intensity   = intensity[intensity > cutOff]
            mz          = np.round(mz, digits)
            mz, intensity= aggregate(mz, intensity)
            yield (mz, intensity)

def readMzXml(path, cutOff=100, digits=2):
    mz, intensity = reduce(merge_runs, get_spectra(path, cutOff, digits))
    return mz, intensity

def readSpectrum(path, cutOff=100, digits=2):
    file_path, ext = os.path.splitext(path)
    reader = {
        'mzxml':readMzXml
    }[ext]
    mz, intensity = reader(cutOff, digits)


# path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
# file_path = path+'FRL_220715_ubi_952_ETD_40ms_04.mzXML'
# file_path = path+'Ubiquitin_ETD_10 ms_1071.mzXML'
# cutOff = 50
# digits = 2
#
# %%time
# mz, intensity = readMzXml(file_path, cutOff, digits )
