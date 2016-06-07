import numpy as np
import os

def openSpectrum(file):
	filename, file_extension = os.path.splitext(file)
	if file_extension in ['.txt']:
		spectrum = []
		for line in open(file):
			l = line.split()
			spectrum.append( tuple( [ float(v) for v in l] ) )
		spectrum = np.array(spectrum)
		return spectrum
	else: 
		raise IOError()

# file = '/Volumes/doom/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodon/MassTodonPy/spectraParsing/Melphalan_UBQ.txt'
# A = openSpectrum(file)
# openSpectrum('dupa.dup')