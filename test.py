%load_ext autoreload
%autoreload
from Formulator.formulator import genMolecules, makeFragments, pandizeSubstances
from Formulator.fasta2atomcnt import fasta2atomCnt
from InSilico.spectrumGenerator import insilicoSpectrum, genIsotopicEnvelope, flatten, makeNoise
from math import sqrt
from Formulator.isotopeCalculator import IsotopeCalculations
from Solver.solver import get_linprog_input
from scipy.optimize import linprog

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
cheb = .01
jointProb 		= .999
precisionDigits = 2
precisionMass  	= .05 # In Daltons; the radius of mass bucket around theoretical peaks

G = nx.Graph()
tolInt 	= it.IntervalTree()
mols 	= genMolecules(fasta, 3, 'cz',modifications)

# Nodes of molecules.
for i, molInfo in enumerate(mols):
	molType, q, g, atomCnt, monoM, meanM, stDevM = molInfo
	mol = ('M',i)
	tolInt[ meanM-stDevM/sqrt(cheb) : meanM+stDevM/sqrt(cheb) ] = mol
	G.add_node( mol, q=q, g=g, atomCnt=atomCnt, massMono=monoM, massMean=meanM, massStDev=stDevM )

# Nodes of experimental peaks; edges beteen experimental peaks and molecules
for i, (mz,intensity) in enumerate(MassSpectrum):
	exp_peak = ('E',i)
	G.add_node( exp_peak, mz=mz, intensity=intensity)
	for tolData in tolInt[mz]:
		molecule = tolData.data
		G.add_edge( molecule, exp_peak, M2E=True )

IC = IsotopeCalculations()
tolInt = it.IntervalTree()
isoCnt = 0

# Isotopologue nodes added
for mol, data in G.nodes_iter(data=True):
	nodeType, nodeNo = mol
	if nodeType=='M' and G.neighbors(mol) > 0:
		for mz, prob in IC.getIsotopicEnvelope( data['atomCnt'], jointProb, precisionDigits):
			iso = ('I',isoCnt)
			G.add_node( iso, mz=mz, prob=prob )
			G.add_edge( mol, iso )
			isoCnt += 1
			tolInt[ mz-precisionMass : mz+precisionMass ] = iso

# Edges between isotopologue nodes and experimental peaks added
for exp_peak, data in G.nodes_iter(data=True):
	nodeType, nodeNo = exp_peak
	if nodeType=='E':
		for tolData in tolInt[data['mz']]:
			iso = tolData.data
			G.add_edge( iso, exp_peak, I2E=True )

# Removing edges between experimental peaks and molecules
G.remove_edges_from( (v1, v2) for v1, v2, d in G.edges_iter(data=True) if 'M2E' in d )

# Form groups of experimental peaks
Gcnt = 0
eGs = dict()

for exp_peak, exp_peakData in G.nodes_iter(data=True):
	nodeType, exp_peakNo = exp_peak
	if nodeType=='E':
		I_neigh = frozenset(G.neighbors(exp_peak))
		if len(I_neigh)>0:
			if not I_neigh in eGs:
				eG = ('eG', Gcnt) # Experiment Group
				eGs[ I_neigh ] = eG
				Gcnt += 1
				G.add_node( eG, intensity = 0.0 )
				for iso in I_neigh:
					G.add_edge( iso, eG )
			eG = eGs[ I_neigh ]
			G.node[ eG ]['intensity'] += exp_peakData['intensity']
			G.add_edge( eG, exp_peak )

# Remove edges between experimental peaks and isotopologues: bridging through groups of experimental peaks only
G.remove_edges_from( (v1, v2) for v1, v2, d in G.edges_iter(data=True) if 'I2E' in d )

def isDeconvoProb(g):
	isProblem = False
	for (n1,id1), (n2,id2) in g.edges_iter():
		if n1 == 'eG' or n2 == 'eG':
			isProblem = True
			break
	return isProblem

from itertools import ifilter

deconvoProbs = ifilter( isDeconvoProb, nx.connected_component_subgraphs(G) )



# dProbs = list(deconvoProbs)
# g = dProbs[3]


X = [ get_linprog_input(g) for g in deconvoProbs]
g, coord2edges, linprogInput = X[0]
linprog(options={"disp": True}, **linprogInput)



res = linprog(c, A_ub, b_ub, A_eq, b_eq, (0.0, None), 'simplex', None, {"disp": True})

res

# X[0][0]
# import matplotlib.pyplot as plt
# nx.draw(X[2][0])
# plt.show()
