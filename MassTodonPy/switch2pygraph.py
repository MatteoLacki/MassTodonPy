%load_ext autoreload
%autoreload
%load_ext line_profiler

from MassTodon import MassTodon

fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 9; modifications = {}; ionsNo  = 10000; P = .999

massTodon = MassTodon(fasta, Q, massPrecDigits = 1)
masses, intensities, noise_masses, noise_intensities = massTodon.randomSpectrum(ionsNo)

mz      = np.append(masses, noise_masses)
ints    = np.append(intensities, noise_intensities)

massSpectrum = list(zip(mz,ints))
massTodon.peakPicker.setMassSpectrum(massSpectrum)

import igraph as ig
import  intervaltree
from math import sqrt


#####

from IsotopeCalculator.isotopeCalculator import atomCnt2string



%%time
tolInt  = intervaltree.IntervalTree()
molecules = []
edges   = []

for mol_cnt, (mol, q, g) in enumerate( massTodon.formulator.makeMolecules()):
    atomCnt = atomCnt2string(mol['atomCnt'])
    name    = 'M'+str(mol_cnt)
    meanMass= massTodon.isoCalc.getMassMean(mol['atomCnt'])
    massVar = massTodon.isoCalc.getMassVar(mol['atomCnt'])
    chebDev = sqrt(massVar/massTodon.peakPicker.cheb)
    tolInt[ meanMass-chebDev : meanMass+chebDev ] = name
    molecules.append({'name':name, 'q':q, 'g':g, 'atomCnt':atomCnt, 'type':'M'})

for E_cnt, (mz,intensity) in enumerate(massTodon.peakPicker.MS):
    E = 'E'+str(E_cnt)
    mol = {'name':name, 'mz':mz, 'intensity':intensity, 'type':'E' }
    for tolData in tolInt[mz]:
        M = tolData.data
        edges.append({ 'source': E, 'target': M})

G = ig.Graph.DictList(vertices=molecules, edges=edges, iterative=False)


# Add isotopologue nodes (I) to the problem graph G.
tolInt  = intervaltree.IntervalTree()
isoCnt  = 0




fP = formulaParser()
fP.parse('C100H200')




f
# _neigborhood_size_gt
for f in G.vs(type='M', _degree_gt=1 ):
    print f

G.neighborhood_size(f)


for mol, data in G.nodes_iter(data=True):
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



#####


# atomCnt = mol['atomCnt']
# meanM   = massTodon.isoCalc.getMassMean(atomCnt)            # Change massTodon to self as it will be used by the peakPicker class
# monoM   = massTodon.isoCalc.getMonoisotopicMass(atomCnt)
# M       = ('M', mol_cnt)
#
# chebDev = sqrt(massVar/massTodon.peakPicker.cheb) # Change to self.
#
# self.G.add_node( M, q=q, g=g, atomCnt=atomCnt, massMono=monoM, massMean=meanM, massStDev=sqrt(massVar) )








massTodon.peakPicker.add_M_n_E()



#TODO what goes wrong with graph generation: why so slow?
massTodon.peakPicker.add_I()
massTodon.peakPicker.add_eG()
G = massTodon.peakPicker.G
