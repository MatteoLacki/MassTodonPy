import intervaltree
import networkx as nx
from Formulator.isotopeCalculator import IsotopeCalculations
from InSilico.spectrumGenerator import genIsotopicEnvelope
from math import sqrt

class peakPicker():
    '''Class for peak picking.'''

    def __init__(self,
            molecules_iter,
            MassSpectrum,
            chebyshevCoverage       = 0.99,
            jointProbabilityIsoSpec = .999,
            precisionDigits         = 2,
            precisionMass           = .05
    ):
        self.cheb   = 1.0 - chebyshevCoverage
        self.G      = nx.Graph()
        self.mols   = molecules_iter
        self.jP     = jointProbabilityIsoSpec
        self.MS     = MassSpectrum
        self.prec   = precisionDigits
        self.massPrec = precisionMass

    def add_M_n_E(self):
        '''Add molecule (M) and experimental (E) peak nodes to the problem graph G and links them using the idea of tolerance interval.'''

        # Nodes of molecules.
        tolInt  = intervaltree.IntervalTree()
        for i, molInfo in enumerate(self.mols):
            molType, q, g, atomCnt, monoM, meanM, stDevM = molInfo
            mol = ('M',i)
            chebDev = stDevM/sqrt(self.cheb)
            tolInt[ meanM-chebDev : meanM+chebDev ] = mol
            self.G.add_node( mol, q=q, g=g, atomCnt=atomCnt, massMono=monoM, massMean=meanM, massStDev=stDevM )

        # Nodes of experimental peaks; edges beteen experimental peaks and molecules
        for i, (mz,intensity) in enumerate(self.MS):
            exp_peak = ('E',i)
            self.G.add_node( exp_peak, mz=mz, intensity=intensity)
            for tolData in tolInt[mz]:
                molecule = tolData.data
                self.G.add_edge( molecule, exp_peak, M2E=True )

    def add_I(self):
        '''Add isotopologue nodes (I) to the problem graph G.'''
        IC      = IsotopeCalculations()
        tolInt  = intervaltree.IntervalTree()
        isoCnt  = 0

        # Add isotopologue nodes I
        for mol, data in self.G.nodes_iter(data=True):
        	nodeType, nodeNo = mol
        	if nodeType=='M' and self.G.neighbors(mol) > 0:
        		for mz, prob in IC.isoEnvelope( data['atomCnt'], self.jP, self.prec):
        			iso = ('I',isoCnt)
        			self.G.add_node( iso, mz=mz, prob=prob )
        			self.G.add_edge( mol, iso )
        			isoCnt += 1
        			tolInt[ mz - self.massPrec : mz + self.massPrec ] = iso

        # Edges between isotopologue nodes and experimental peaks added
        for exp_peak, data in self.G.nodes_iter(data=True):
        	nodeType, nodeNo = exp_peak
        	if nodeType=='E':
        		for tolData in tolInt[data['mz']]:
        			iso = tolData.data
        			self.G.add_edge( iso, exp_peak, I2E=True )

        # Removing edges between experimental peaks and molecules
        self.G.remove_edges_from(
            (v1, v2) for v1, v2, d in self.G.edges_iter(data=True) if 'M2E' in d )

    def add_eG(self):
        '''Add experimental group nodes (eG) to the problem graph G.'''
        Gcnt= 0
        eGs = dict()
        for E, E_data in self.G.nodes_iter(data=True):
        	nodeType,_ = E
        	if nodeType=='E':
        		I_neigh = frozenset(self.G.neighbors(E))
        		if len(I_neigh)>0:
        			if not I_neigh in eGs:
        				eG = ('eG', Gcnt) # Experiment Group
        				eGs[ I_neigh ] = eG
        				Gcnt += 1
        				self.G.add_node( eG, intensity = 0.0 )
        				for iso in I_neigh:
        					self.G.add_edge( iso, eG )
        			eG = eGs[ I_neigh ]
        			self.G.node[ eG ]['intensity'] += E_data['intensity']
        			self.G.add_edge( eG, E )

        # Remove edges between experimental peaks and isotopologues: bridging through groups of experimental peaks only
        self.G.remove_edges_from(
            (v1, v2) for v1, v2, d in self.G.edges_iter(data=True) if 'I2E' in d )

    def pickPeaks(self):
        '''Perform peak picking.'''
        self.add_M_n_E()
        self.add_I()
        self.add_eG()
        return self.G
