from spectrumPreprocessing.readMzXML import openMzXML 

spectra, metaInfo = openMzXML("~/Documents/Science/MassTodon/MassTodonPy/data/FRL_220715_ubi_952_ETD_40ms_01.mzXML")
threshold_on_intensity = 10.

spectra 	= [ [ (m, i) for m, i in s if i > threshold_on_intensity ] for s in spectra ]
spectrum 	= spectra[0]


