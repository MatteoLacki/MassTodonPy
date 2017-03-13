%load_ext autoreload
%autoreload
%load_ext line_profiler

from MassTodon import MassTodon
import numpy as np
import igraph as ig
from intervaltree import Interval as I, IntervalTree as Itree
from math import sqrt
from pandas import DataFrame
from frozendict import frozendict
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
from data.mzXMLr2py import openMzXML


fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; Prob = .999; mzPrec = .1

res = openMzXML(file='/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/FRL_220715_ubi_952_ETD_40ms_04.mzXML')



massTodon = MassTodon(fasta, Q, massPrecDigits = 1)

masses, intensities = massTodon.isoCalc.randomPrecursors(fasta,Q,ionsNo,massTodon.formulator)
noise_masses, noise_intensities = massTodon.isoCalc.addNoise( masses, intensities, .05 )

mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)

massSpectrum = list(zip(mz,ints))
massTodon.peakPicker.setMassSpectrum(massSpectrum)

DataFrame(massSpectrum).to_csv(path_or_buf='/Users/matteo/Documents/MassTodon/insilicoSpectrum.csv', index=False)

# len(massSpectrum)


### approach focussed on Divide ed Impera
ePeaks = Itree( I(mz-mzPrec, mz+mzPrec, (mz, intensity)) for mz, intensity in massTodon.peakPicker.MS )
G      = nx.Graph()

for cnt, (mType, formula, aaNo, q, g) in enumerate( massTodon.formulator.makeMolecules() ):
    massMono, massMean, massVar = massTodon.isoCalc.getSummary(formula)
    mzMean  = (float(massMean)+g+q)/q
    chebDev = sqrt(massVar/massTodon.peakPicker.cheb)
    eAround = ePeaks[mzMean]
    if eAround:
        M = 'M'+str(cnt)
        G.add_node( M, q=q, g=g, formula=formula )
        for e in eAround:
            Emz, Eintensity = e.data
            G.add_node( Emz, intensity=Eintensity )
            G.add_edge( M, Emz, M2E=True )

ConnectedComponents = list(nx.connected_component_subgraphs(G))
len(ConnectedComponents)

Counter([ len(cc) for cc in ConnectedComponents ])

for cc in ConnectedComponents:
    if len(cc) == 16:
        X = cc

# nx.draw(X)
# plt.show()








X = ConnectedComponents[1110]
nx.draw(X)
plt.show()

### approach focussed on intervaltrees - wrong
%%time
mzPrec = .05
ePeaks = Itree( I(mz-mzPrec, mz+mzPrec, intensity) for mz, intensity in massTodon.peakPicker.MS )

iPeaks  = Itree()
IsoPeaks= {}
edges   = []
vertices= []


iCnt = 0
for cnt, (mType, formula, aaNo, q, g) in enumerate( massTodon.formulator.makeMolecules() ):
    massMono, massMean, massVar = massTodon.isoCalc.getSummary(formula)
    mzMean  = (float(massMean)+g+q)/q
    chebDev = sqrt(massVar/massTodon.peakPicker.cheb)
    eAround = ePeaks[mzMean]
    if eAround:
        M = 'M'+str(cnt)
        masses, probs = massTodon.isoCalc.isoEnvelope( formula, Prob, q, g )
        for mz_iso, pr_iso in zip(masses,probs):
            ePeaks[mz_iso]




masses, probs


test = Itree()
test.update(I(i,i+1,i**2) for i in xrange(5))

ePeaks[31.3]
ePeaks
%%time
x = np.dstack((masses,probs))

%%time
x = np.array(zip(masses,probs))


### Old approach
%%time
tolInt      = intervaltree.IntervalTree()
molecules   = []
edges       = []

for mol_cnt, (molType, atomCnt_str, sideChainsNo, q, g) in enumerate( massTodon.formulator.makeMolecules()):
    M = 'M'+str(mol_cnt)
    monoMass, meanMass, massVar = massTodon.isoCalc.getSummary(atomCnt_str)
    meanMass= (float(meanMass)+g+q)/q
    chebDev = sqrt(massVar/massTodon.peakPicker.cheb)
    tolInt[ meanMass-chebDev : meanMass+chebDev ] = M
    # G.add_vertex( type='M', name=M, q=q, g=g, atomCnt_str=atomCnt_str ) # slower
    molecules.append({'name':M, 'q':q, 'g':g, 'atomCnt_str':atomCnt_str, 'type':'M'})

for E_cnt, (mz,intensity) in enumerate(massTodon.peakPicker.MS):
    E = 'E'+str(E_cnt)
    M_around_E = tolInt[mz]
    if M_around_E:
        molecules.append({'name':E, 'type':'E', 'mz':mz, 'intensity':intensity})
        # G.add_vertex( type='E', name=E, mz=mz, intensity=intensity ) # slower
        for tolData in tolInt[mz]:
            M = tolData.data
            edges.append({'source':E, 'target':M, 'M2E':True })

G = ig.Graph.DictList(vertices=molecules, edges=edges, iterative=False)
# G = ig.Graph.DictList(vertices=molecules, edges=edges, iterative=False)

list(G.vs(type='E'))
list(G.es())
# Add isotopologue nodes (I) to the problem graph G.
tolInt  = intervaltree.IntervalTree()
isoCnt  = 0

## Testing interval trees capability to cope with experimental peaks.
massPrec = .05

for E_cnt, (mz,intensity) in enumerate(massTodon.peakPicker.MS):

%%time
experimental_tree = IntervalTree( Interval(mz-massPrec, mz+massPrec, intensity) for mz, intensity in massTodon.peakPicker.MS ) # Not so fast..




# For every molecule paired with any experimental peaks (can enhance to peaks with height above some value)
#   check if the experimental peaks around
for mol, data in self.G.nodes_iter(data=True):
    nodeType, nodeNo = mol
    if (nodeType=='M') and (self.G.neighbors(mol) > 0):
        atomCnt = data['atomCnt']
        q, g = data['q'], data['g']
        masses, probs = self.isoCalc.isoEnvelope(atomCnt, self.jP, q, g)
        for mz, prob in zip(masses, probs):
            iso = ('I',isoCnt)
            self.G.add_node( iso, mz=mz, prob=prob )
            self.G.add_edge( mol, iso )
            isoCnt += 1
        tolInt[ mz - self.massPrec : mz + self.massPrec ] = iso
