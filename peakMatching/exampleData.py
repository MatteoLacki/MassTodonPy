savePath = '/Users/matteo/Dropbox/Science/MassSpectrometry/MassTodonPy'
tol = .05
minExpSupp = .6

theory = {
	'theory_0' : [(20.1, 1.), (21.1, 2.), (22.2, 3.)],
	'theory_1' : [(21.,.1),(22.,.2),(23.,.3),(24.,5.),(25.,10.),(26.,6.)],
	'theory_2' : [(22.1,.1),(23.11,.5),(24.14,.8),(25.13,1.)],
	'theory_3' : [(40.,.2),(41.,.3),(42.,.1)],
	'theory_4' : [(41.01,.1),(42.01,.6),(43.01,.7),(44.01,.01)]
}
# empiria = [
# 	(20.1, 10.), (21.101, 20.),(22.09, 32.),(24.15, 10),(25.1, 14),
# 	(40.02,213.),(41.02,521.),(42.01,122.),
# 	(41.09,100.),(42.,640.),(43.02, 723),
# 	(41.051,10.),
# 	(50,1.)
# ]

def normalise(envelope):
	'Normalising theoretical weights to probabilities'
	totWeight = sum( float(weight) for mass, weight in envelope )
	return [ (mass, weight/totWeight) for mass, weight in envelope ]

for compound in theory:
	theory[compound] = normalise( theory[compound] )

ionCounts={ 'theory_0': 100.,
			'theory_1': 200.,
			'theory_2': 300.,
			'theory_3': 400.,
			'theory_4': 500. 	}	

empiria = []
for th in theory:
	ioNo = ionCounts[th]
	for mass, intensity in theory[th]:
		empiria.append(( mass, intensity*ioNo ))