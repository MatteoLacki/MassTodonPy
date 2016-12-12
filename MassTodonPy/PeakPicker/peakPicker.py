# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.

import  intervaltree
import  networkx        as      nx
from    math            import  sqrt

class PeakPicker():
    '''Class for peak picking.'''

    def __init__(self,
            formulator,
            isotopeCalculator,
            chebyshevCoverage       = 0.99,
            jointProbabilityIsoSpec = .999,
            precisionDigits         = 2,
            precisionMass           = .05   ):

        self.cheb   = 1.0 - chebyshevCoverage
        self.G      = nx.Graph()

        self.formulator = formulator
        self.isoCalc = isotopeCalculator

        self.jP     = jointProbabilityIsoSpec
        self.prec   = precisionDigits
        self.massPrec = precisionMass

    def setMassSpectrum(self, massSpectrum):
        self.MS = massSpectrum

    def add_M_n_E(self):
        '''Add molecule (M) and experimental (E) peak nodes to the problem graph G and links them using the idea of tolerance interval.'''

        # Nodes of molecules.
        tolInt  = intervaltree.IntervalTree()
        for mol_cnt, (mol, q, g) in enumerate(self.formulator.makeMolecules()):
            atomCnt = mol['atomCnt']
            meanM   = self.isoCalc.getMassMean(atomCnt)
            monoM   = self.isoCalc.getMonoisotopicMass(atomCnt)
            M       = ('M', mol_cnt)
            massVar = self.isoCalc.getMassVar(atomCnt)
            chebDev = sqrt(massVar/self.cheb)
            tolInt[ meanM-chebDev : meanM+chebDev ] = M
            self.G.add_node( M, q=q, g=g, atomCnt=atomCnt, massMono=monoM, massMean=meanM, massStDev=sqrt(massVar) )

        # Nodes of experimental peaks; edges beteen experimental peaks and molecules
        for E_cnt, (mz,intensity) in enumerate(self.MS):
            E = ('E', E_cnt) # experimental peak
            self.G.add_node( E, mz=mz, intensity=intensity )
            for tolData in tolInt[mz]:
                M = tolData.data
                self.G.add_edge( M, E, M2E=True )

    def add_I(self):
        '''Add isotopologue nodes (I) to the problem graph G.'''
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
        print 'done'

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
