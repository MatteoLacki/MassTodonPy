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
		if types==alphaEdge or types==xEdge:
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

def solveProblem(g):
	g, coord2edges, linprogInput = get_linprog_input(g)
	res = linprog(**linprogInput)
	return res, g, coord2edges, linprogInput

def solutions_iter(deconvolutionProblems_iter):
	for g in deconvolutionProblems_iter:
		yield solveProblem(g)
