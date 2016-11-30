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

# Nodes of molecules.
for i, molInfo in enumerate(molecules):
	molType, q, g, atomCnt, massMono, massMean, massStDev = molInfo
	tag = ('molecule',i)
	tolIntervals[ massMean-massStDev/sqrt(chebyshevEps) : massMean+massStDev/sqrt(chebyshevEps) ] = tag
	G.add_node( tag, q=q, g=g, atomCnt=atomCnt, massMono=massMono, massMean=massMean, massStDev=massStDev )

# Nodes of experimental peaks
# Edges beteen experimental peaks and molecules
for i, (mass,intensity) in enumerate(MassSpectrum):
	experimentalPeak = ('experimental',i)
	G.add_node( experimentalPeak, mz=mass, I=intensity)
	for tolData in tolIntervals[mass]:
		molecule = tolData.data
		G.add_edge( experimentalPeak, molecule, experiment2molecule=True)

IC = IsotopeCalculations()
tolIntervals 	= it.IntervalTree()
isotopologueCnt = 0

# Isotopologue nodes added
for node, data in G.nodes_iter(data=True):
	nodeType, nodeNo = node
	if nodeType=='molecule' and G.neighbors(node) > 0:
		for mz, prob in IC.getIsotopicEnvelope(data['atomCnt'], jointProb, precisionDigits):
			isotopologue = ('isotopologue',isotopologueCnt)
			G.add_node( isotopologue, mz=mz, prob=prob)
			G.add_edge( isotopologue, node)
			isotopologueCnt += 1
			tolIntervals[ mz-precisionMass : mz+precisionMass ] = isotopologue

# Edges between isotopologue nodes and experimental peaks added
for node, data in G.nodes_iter(data=True):
	nodeType, nodeNo = node
	if nodeType=='experimental':
		for tolData in tolIntervals[data['mz']]:
			isotopologueTag = tolData.data
			G.add_edge(node, isotopologueTag, experiment2isotopologue=True)

# Removing edges between experimental peaks and molecules
G.remove_edges_from( (v1, v2) for v1, v2, d in G.edges_iter(data=True) if 'experiment2molecule' in d )

# Form groups of experimental peaks
Gcnt = 0
experimentalGroups = dict()
for experimentalPeak, experimentalPeakData in G.nodes_iter(data=True):
	nodeType, experimentalPeakNo = experimentalPeak
	if nodeType=='experimental':
		nbrs = frozenset(G.neighbors(experimentalPeak))
		if len(nbrs)>0:
			if not nbrs in experimentalGroups:
				experimentalGroup = ('experimentalGroup', Gcnt)
				experimentalGroups[nbrs]= experimentalGroup
				Gcnt += 1
				G.add_node( experimentalGroup, I=0.0 )
			experimentalGroup = experimentalGroups[nbrs]
			G.node[experimentalGroup]['I'] += experimentalPeakData['I']
			G.add_edge( experimentalPeak, experimentalGroup )

# Remove edges between experimental peaks and isotopologues: bridging through groups of experimental peaks only
G.remove_edges_from( (v1, v2) for v1, v2, d in G.edges_iter(data=True) if 'experiment2isotopologue' in d )

len(list(nx.connected_component_subgraphs(G)))

len(G.nodes(data=True))
len(G.edges(data=True))

experimentalGroups


import matplotlib.pyplot as plt
nx.draw(G)
plt.show()
G.edges(data=True)

tolIntervals
MassSpectrum

molecules
