%load_ext autoreload
%autoreload
%load_ext line_profiler

from MassTodon import MassTodon
import numpy as np
import igraph as ig
import  intervaltree
from math import sqrt
from pandas import DataFrame
from frozendict import frozendict
from collections import Counter

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; P = .999

massTodon = MassTodon(fasta, Q, massPrecDigits = 1)
masses, intensities, noise_masses, noise_intensities = massTodon.randomSpectrum(ionsNo)
mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)
massSpectrum = list(zip(mz,ints))
massTodon.peakPicker.setMassSpectrum(massSpectrum)

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
            edges.append({'source':E, 'target':M})

G = ig.Graph.DictList(vertices=molecules, edges=edges, iterative=False)







list(G.vs(type='E'))
# Add isotopologue nodes (I) to the problem graph G.
tolInt  = intervaltree.IntervalTree()
isoCnt  = 0


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
