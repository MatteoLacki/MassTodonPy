import intervaltree as iTree
 
tree 	= iTree.IntervalTree()
tol 	= .05
theory 	= [
	[(20.1, 1.), (21.1, 2.), (22.2, 3.)],
	[(21.,.1),(22.,.2),(23.,.3),(24.,5.),(25.,10.),(26.,6.)],
	[(22.1,.1),(23.11,.5),(24.14,.8),(25.13,1.)],
	[(40.,.2),(41.,.3),(42.,.1)],
	[(41.01,.1),(42.01,.6),(43.01,.7)]
]

empiria = [
	(20.1, 10.), (21.101, 20.),(22.09, 32.),(24.15, 10),(25.1, 14),
	(40.02,213.),(41.02,521.),(42.01,122.),
	(41.09,100.),(42.,640.),(43.02, 723)
]

for famNo, spectrum in enumerate(theory):
	for peakNo, peak in enumerate(spectrum):
		mass, intensity = peak
		tree[ mass-tol : mass+tol ] = (peakNo, famNo)


for mass, intensity in empiria:
	print tree[mass]
	for x in tree[mass]:
		print x.data