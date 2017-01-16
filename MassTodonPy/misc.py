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

from collections import defaultdict
from linearCounter import linearCounter as lCnt
from itertools import ifilter
from networkx import connected_component_subgraphs as connectedComponents


def standardize(modifications):
    '''Standardize modifications so that they meet the internal nomenclature scheme.'''
    backboneAtom2aaNomen = {'N':'L', 'Calpha':'C', 'C':'R'}
    R = defaultdict(lambda:defaultdict(lCnt))
    for tag, atomCnt in modifications.items():
        R[ tag[1]-1 ][ backboneAtom2aaNomen[tag[0]] ] = lCnt(atomCnt)
    return R

def countIsNegative(atomCnt):
    '''Check if any element of a dictionary is a negative number.'''
    return any( atomCnt[elem]<0 for elem in atomCnt )

def toDeconvolve(g):
	isProblem = False
	for (n1,id1), (n2,id2) in g.edges_iter():
		if n1 == 'eG' or n2 == 'eG':
			isProblem = True
			break
	return isProblem

def deconvIter(G):
    return ifilter(toDeconvolve, connectedComponents(G))

def atomCnt2string(atomCnt):
    '''Translate a dictionary of atom counts into a uniquely defined string.'''
    keys = atomCnt.keys()
    keys.sort()
    return "".join( el+str(atomCnt[el]) for el in keys )


# class Indexer(dict):
# 	def get(self, v1, v2):
# 		args = frozenset( (v1,v2) )
# 		return super(Indexer, self).__getitem__(args)
# 	def add(self, v1, v2, val):
# 		args = frozenset( (v1,v2) )
# 		super(Indexer, self).__setitem__(args, val)
