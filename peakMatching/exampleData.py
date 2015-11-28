savePath = '/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodonPy'
tol = .05
minExpSupp = .6

theory 	= [
	[(20.1, 1.), (21.1, 2.), (22.2, 3.)],
	[(21.,.1),(22.,.2),(23.,.3),(24.,5.),(25.,10.),(26.,6.)],
	[(22.1,.1),(23.11,.5),(24.14,.8),(25.13,1.)],
	[(40.,.2),(41.,.3),(42.,.1)],
	[(41.01,.1),(42.01,.6),(43.01,.7),(44.01,.01)]
]
# empiria = [
# 	(20.1, 10.), (21.101, 20.),(22.09, 32.),(24.15, 10),(25.1, 14),
# 	(40.02,213.),(41.02,521.),(42.01,122.),
# 	(41.09,100.),(42.,640.),(43.02, 723),
# 	(41.051,10.),
# 	(50,1.)
# ]

def normalise(envelope):
	totWeight = sum( float(weight) for mass, weight in envelope )
	return [ (mass, weight/totWeight) for mass, weight in envelope ]
	# Normalising weight to probabilities
theory = [ normalise(spectrum) for spectrum in theory ]

#Check if it is not better to include these peaks using IntervalTree( <iterable> )
	# IntervalTree( Interval( mass-tol, mass+tol, (peakNo,famNo) ) for  )
	
ionNo = [ 100., 200., 300., 500., 1000. ]

empiria = [
	(mass, intensity * ioNo) for spectrum, ioNo in zip(theory, ionNo) for mass, intensity in spectrum
]
empiria.append(( 100.,200. ))
empiria.append(( 20.11, 200. ))
# print empiria

