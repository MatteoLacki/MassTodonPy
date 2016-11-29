%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules, makeFragments, pandizeSubstances
from Formulator.fasta2atomcnt import fasta2atomCnt
from InSilico.spectrumGenerator import insilicoSpectrum, genIsotopicEnvelope, flatten, makeNoise
from math import sqrt
from collections import Counter
from Formulator.isotopeCalculator import IsotopeCalculations

fasta = substanceP  = 'RPKPQQFFGLM'
Q = 3; modifications = {}; ionsNo = 1000000; P = .999

IS = insilicoSpectrum(fasta, Q, ionsNo, P)
sample, probs 	= IS.rvs(ionsNo)
MassSpectrum 	= [ (m, float(i))for m, i in flatten(sample)]

Noise = makeNoise( MassSpectrum, percentPeaks = .2 )
MassSpectrum.extend(Noise)
########################################################################
from operator import itemgetter
MassSpectrum.sort(key=itemgetter(0))

import igraph as ig
import intervaltree as it
import networkx as nx
########################################################################
chebyshevEps 	= .01
jointProb 		= .999
precisionDigits = 3
precisionMass  	= .05 # In Daltons; the radius of mass bucket around theoretical peaks

G = nx.Graph()
tolIntervals= it.IntervalTree()
molecules 	= genMolecules(fasta, 3, 'cz',modifications)

for i, molInfo in enumerate(molecules):
	molType, q, g, atomCnt, massMono, massMean, massStDev = molInfo
	tag = ('molecule', str(i))
	tolIntervals[ massMean-massStDev/sqrt(chebyshevEps) : massMean+massStDev/sqrt(chebyshevEps) ] = tag
	G.add_node( tag, q=q, g=g, atomCnt=atomCnt, massMono=massMono, massMean=massMean, massStDev=massStDev )

for i, (mass,intensity) in enumerate(MassSpectrum):
	tag = ('experimental',str(i))
	G.add_node(tag)
	for tolData in tolIntervals[mass]:
		moleculeTag = tolData.data
		G.add_edge(tag, moleculeTag, filtering=0)

IC = IsotopeCalculations()

tolIntervals 	= it.IntervalTree()
isotopologueCnt = 0
for node, data in G.nodes_iter(data=True):
	nodeType, nodeNo = node
	if nodeType=='molecule' and G.neighbors(node) > 0:
		for mz, prob in IC.getIsotopicEnvelope(data['atomCnt'], jointProb, precisionDigits):
			tag = ('isotopologue',isotopologueCnt)
			G.add_node(tag, mz=mz, prob=prob)
			G.add_edge(tag, node, filtering=1)
			isotopologueCnt += 1
			tolIntervals[ mz-precisionMass : mz+precisionMass ] = tag

for node, data in G.nodes_iter(data=True):
	nodeType, nodeNo = node
	if nodeType=='experimental' and G.neighbors(node) > 0:
		for tolData in tolIntervals[mass]:
			isotopologueTag = tolData.data
			G.add_edge(node, isotopologueTag, filtering=1)

G.edges(data=True)


G.edges()

len(G.edges())



import matplotlib.pyplot as plt
nx.draw(G)
plt.show()
G.edges(data=True)

tolIntervals
MassSpectrum

molecules
