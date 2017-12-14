%load_ext autoreload
%autoreload 2

import matplotlib.pyplot as plt
import networkx as nx
from networkx import connected_component_subgraphs as connected_components
from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Data.Constants import eps
from MassTodonPy.Deconvolutor.Deconvolutor import deconvolve
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch, diag, incidence_matrix
from MassTodonPy.MatchMaker.CzMatch import CzMatch
from MassTodonPy.MatchMaker.RegularizedCzMatch import RegularizedCzMatch

# %%time

# twice bigger substanceP molecule.

mol = get_dataset('substanceP') # adjust the spectrum
mol.precursor.q = 4
mol.precursor.fasta += mol.precursor.fasta
mol.precursor.formula += mol.precursor.formula
mols = list(mol.precursor.molecules())

mols_dict = {(m.name, m.q, m.g): m for m in mols}
results = [
    {'alphas': [
        {'estimate': 1000,
         'molecule': mols_dict[('c10', 1, 0)]},
        {'estimate': 1300,
         'molecule': mols_dict[('c10', 2, 0)]},
        {'estimate': 5000,
         'molecule': mols_dict[('c5', 1, 0)]},
        {'estimate': 8000,
         'molecule': mols_dict[('c6', 1, 0)]}],
     'status': 'optimal'},
    {'alphas': [
        {'estimate': 1200,
         'molecule': mols_dict[('z12', 1, 0)]},
        {'estimate': 4800,
         'molecule': mols_dict[('z6', 1, 0)]},
        {'estimate': 8300,
         'molecule': mols_dict[('z5', 1, 0)]}],
     'status': 'optimal'},
    {'alphas': [
        {'estimate': 40000,
         'molecule': mols_dict[('precursor', 4, 0)]},
        {'estimate': 40000,
         'molecule': mols_dict[('precursor', 3, 0)]},
        {'estimate': 30000,
         'molecule': mols_dict[('precursor', 2, 0)]},
        {'estimate': 20000,
         'molecule': mols_dict[('precursor', 1, 0)]}],
     'status': 'optimal'}]

simple_matches = SimpleCzMatch(results, mol.precursor)

simple_matches.get_intensities()
simple_matches.get_probabilities()
simple_matches.get_branching_ratios()

matches = CzMatch(results, mol.precursor)
matches.get_intensities()
matches.get_probabilities()
matches.get_branching_ratios()

Q = matches.precursor.q
G = matches.graph

import numpy as np
# lavish: all fragments lose cofragments
lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
node_intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
costs = matrix([float(ETnoD_PTR) for N,M, ETnoD_PTR in G.edges.data('ETnoD_PTR')])
edges_cnt = G.size()
equalities = incidence_matrix(G, len(node_intensities), edges_cnt)
inequalities = diag(-1.0, edges_cnt)
upper_bounds = matrix([0.0] * edges_cnt)

print(costs, equalities, node_intensities, inequalities, upper_bounds)
primalstart = {}
primalstart['x'] = matrix([0.0] * edges_cnt)
primalstart['s'] = matrix([eps] * len(upper_bounds))

solution = solvers.conelp(c=costs,
                          G=inequalities,
                          h=upper_bounds,
                          A=equalities,
                          b=node_intensities,
                          primalstart=primalstart)
solution['primal objective']

equalities

# removing self-loops
G.remove_edges_from([(N, M) for N, M in G.edges if N is M])
G.remove_nodes_from([N for N, deg in G.degree if not deg])

G.degree

for N, deg in G.degree():
    print(N, deg)
lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
node_intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
PTR_costs = matrix([float(PTR) for N,M, PTR in G.edges.data('PTR')])
edges_cnt = G.size()


inequalities = sparse([incidence_matrix(G, len(node_intensities), edges_cnt),
                       diag(-1.0, edges_cnt)])
upper_bounds = matrix([0.0] * edges_cnt)
inequality_constraints = matrix([node_intensities, upper_bounds])







print(node_intensities)








costs_no_self_loops = matrix([float(PTR) for N, M, PTR in
                              G.edges.data('PTR') if not N is M])
print(costs_no_self_loops)
node_intensities_no_self_loops = matrix([float(I) for N, I
                                         in G.nodes.data('intensity')
                                         if G.degree[N] > 2])

node_intensities_no_self_loops

MR_equalities = sparse([equalities, costs.T])
print(equalities, MR_equalities)
ETnoD_costs = matrix([float(ETnoD) for N,M, ETnoD_PTR in G.edges.data('ETnoD_PTR')])










for i, (N, M) in enumerate(G.edges()):
    G[N][M]['flow'] = solution['x'][i]





G.edges.data('flow')
x = np.array(solution['x'])
np.array(solution['x'][0:4]).sum()

from itertools import islice



for N, M, ETnoD_PTR in islice(G.edges.data('ETnoD_PTR'), 4):
    print(N, M, ETnoD_PTR, G[N][M]['flow'])
sum(ETnoD_PTR * G[N][M]['flow'] for N, M, ETnoD_PTR in islice(G.edges.data('ETnoD_PTR'), 4))


print(x[0:4])
print(solution['x'][0:4].T * matrix([1,2,0,1]))

G.edges.data('ETnoD_PTR')
print(np.round(solution['x']))
sum(I)





# matches.get_probabilities()
# matches.get_intensities()
# matches = CzMatch(results, mol.precursor)
# matches.get_probabilities()
# matches.get_intensities()
#
# matches.graph.edges.data('flow')
