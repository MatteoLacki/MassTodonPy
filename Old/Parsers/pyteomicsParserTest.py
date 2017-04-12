from pyteomics import mzxml, auxiliary

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
file = path+'FRL_220715_ubi_952_ETD_40ms_04.mzXML'
file = '/Users/matteo/Documents/MassTodon/ReadingMzxml/data.mzXML'
file = path+'Ubiquitin_ETD_10 ms_1071.mzXML'


spectra = []
with mzxml.read(file) as reader:
    for spectrum in reader:
        mz = spectrum['m/z array']
        intensity = spectrum['intensity array']
        spectra.append((mz,intensity))

len(spectra)

import rpy2.robjects as r
import numpy as np
r.r.source('getAllSpectra.R', verbose=False)

getRspectra	= r.r['getRspectra']
Rres = getRspectra(file)

len(Rres)           # runs
len(Rres[0])        # spectrum and meta Data
len(Rres[0][0])     # mass and intensity of a spectrum
len(Rres[0][0][0])  # vector of masses

def orderDataFromR(x):
    res = [ (np.array(x[i][0][0]), np.array(x[i][0][1])) for i in xrange(len(x)) ]
    return res
R = orderDataFromR(Rres)

def test_equality(j):
    return sum( np.sum(spectra[i][j]-R[i][j]) for i in xrange(len(spectra)) )


cut_off_intensity=50
round_mass_digits=2

with mzxml.read(file) as reader:
    for spectrum in reader:
        mz = spectrum['m/z array']
        intensity = spectrum['intensity array']
