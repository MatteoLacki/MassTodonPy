from heapq import heappush, heappop
from itertools import chain
from operator import itemgetter
from collections import Counter

empiria = [
	(20.1, 10.), (21.101, 20.),(22.09, 32.),(24.15, 10),(25.1, 14),
	(40.02,213.),(41.02,521.),(42.01,122.),
	(41.09,100.),(42.,640.),(43.02, 723)
]

theory = [
	[(20.1, 1.), (21.1, 2.), (22.2, 3.)],
	[(21.,.1),(22.,.2),(23.,.3),(24.,5.),(25.,10.),(26.,6.)],
	[(22.1,.1),(23.11,.5),(24.14,.8),(25.13,1.)],
	[(40.,.2),(41.,.3),(42.,.1)],
	[(41.01,.1),(42.01,.6),(43.01,.7)]
]

famNo = len(theory)
tol = .05
peakNumbers = [ len(spectrum) for spectrum in theory ]

def normalise(envelope, envelopeNo):
	totWeight = sum( float(weight) for mass, weight in envelope )
	return [ (mass, weight/totWeight, envelopeNo) for mass, weight in envelope ]
	
theory = [ normalise(envelope, envelopeNo) for envelopeNo, envelope in enumerate(theory) ]

theory = list(chain(*theory)) 		# New trick	for me
theory.sort( key=itemgetter(0) ) 	# New trick for me

totProb = [.0]*famNo
theCnt = 0
empCnt = 0
massTheory, probTheory, familyTheory = theory[ theCnt ]
massEmpiria, intensityEmpiria = empiria[ empCnt ]

try:
	while theCnt < len(theory) or empCnt < len(empiria):
		if massEmpiria < massTheory - tol:
			empCnt += 1
			massEmpiria, intensityEmpiria = empiria[ empCnt ]
		else:
			if massEmpiria < massTheory + tol:
				totProb[ familyTheory ] += probTheory
				empCnt += 1
				massEmpiria, intensityEmpiria = empiria[ empCnt ]
			else:			
				theCnt += 1
				massTheory, probTheory, familyTheory = theory[ theCnt ]
except IndexError:
	pass

thresholdProb = .6
remainingFamilies = [ fam for fam, prob in enumerate(totProb) if prob > thresholdProb]

theory = [ (mass, prob, famNo) for mass, prob, famNo in theory if famNo in remainingFamilies ]
print theory

# Now, rerun the soft on theoretical peaks outside the range. 
# These should be now the good families, ordered by the beginnings of the intervals.

theCnt = 0
empCnt = 0
massTheory, probTheory, familyTheory = theory[ theCnt ]
activeEnvelopes = set()
massEmpiria, intensityEmpiria = empiria[ empCnt ]
atheoretic = Counter()


# This procedure is simply wrong....
try:
	while theCnt < len(theory) or empCnt < len(empiria):
		if massEmpiria < massTheory - tol:
			empCnt += 1
			massEmpiria, intensityEmpiria = empiria[ empCnt ]
		else:
			if massEmpiria < massTheory + tol:
				activeEnvelopes.add( famNo )
				empCnt += 1
				massEmpiria, intensityEmpiria = empiria[ empCnt ]
			else:			
				theCnt += 1
				massTheory, probTheory, familyTheory = theory[ theCnt ]
except IndexError:
	pass		



