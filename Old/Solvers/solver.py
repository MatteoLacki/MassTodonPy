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
from scipy.optimize import linprog

def get_linprog_input(g):
	'''Prepares the input for the scipy linprog solver.'''
	R = ('R',0) 	# the root node
	g.add_node(R)  	# root is linked to molecules, M
	for nodeType, nodeNo in g.nodes_iter():
		if nodeType=='M':
			g.add_edge(R, (nodeType,nodeNo))

	edges2coord = {} # edges to solver output
	coord2edges = [] # solver output to edges
	varNo = 0
	c = [] # costs

	alphaEdge 	= frozenset(('R','M'))
	xEdge 		= frozenset(('I','eG'))
	for (t1,n1), (t2,n2) in g.edges_iter():
		types = frozenset((t1, t2))
		if types == alphaEdge or types == xEdge:
			coord2edges.append(((t1,n1), (t2,n2)))
			edges2coord[ frozenset(((t1,n1), (t2,n2))) ] = varNo
			if types==alphaEdge:
				c.append(0.0)
			else:
				c.append(-1.0)
			varNo += 1

	A_eq = []
	A_ub_tmp = defaultdict(lambda: [0]*varNo)
	b_ub_tmp = {}

	for M in g.edge[R]:
		for I in g.edge[M]:
			if len(g.edge[I]) > 1:
				if not I == R:
					A_eq_row = [.0] * varNo
					for eG in g.edge[I]:
						if not eG == M:
							IeGidx = edges2coord[frozenset((I,eG))]
							A_eq_row[ IeGidx ] = 1.0
							A_ub_tmp[eG][ IeGidx ] = 1.0
							b_ub_tmp[eG] = g.node[eG]['intensity']
					A_eq_row[ edges2coord[frozenset((R,M))] ] = -g.node[I]['prob']
					A_eq.append(A_eq_row)

	b_eq = [0]*len(A_eq)
	A_ub = []; b_ub = []

	for eG in b_ub_tmp:
		A_ub.append(A_ub_tmp[eG])
		b_ub.append(b_ub_tmp[eG])

	# coord2edges will be necessary to translate the results back
	linprogInput = {'c':c, 'A_eq':A_eq, 'b_eq':b_eq, 'A_ub':A_ub, 'b_ub':b_ub}
	return g, coord2edges, linprogInput

def isDeconvoProb(g):
	isProblem = False
	for (n1,id1), (n2,id2) in g.edges_iter():
		if n1 == 'eG' or n2 == 'eG':
			isProblem = True
			break
	return isProblem


def solveProblem(g):
	'''Apply the simplex algorithm to solve the modified max flow problem.'''
	g, coord2edges, linprogInput = get_linprog_input(g)
	res = linprog(**linprogInput)
	return res, g, coord2edges, linprogInput


def solutions_iter(deconvolutionProblems_iter):
	'''Iterate over all solutions.'''
	for g in deconvolutionProblems_iter:
		yield solveProblem(g)
