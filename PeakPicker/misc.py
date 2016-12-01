from itertools import ifilter
from networkx import connected_component_subgraphs as connectedComponents

def toDeconvolve(g):
	isProblem = False
	for (n1,id1), (n2,id2) in g.edges_iter():
		if n1 == 'eG' or n2 == 'eG':
			isProblem = True
			break
	return isProblem

def deconvIter(G):
    return ifilter(toDeconvolve, connectedComponents(G))
