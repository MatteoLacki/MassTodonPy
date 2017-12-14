%load_ext autoreload
%autoreload 2

from cvxopt import  matrix, spmatrix, sparse, spdiag, solvers
import matplotlib.pyplot as plt
import networkx as nx
from networkx import connected_component_subgraphs as connected_components
import numpy as np

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Data.Constants import eps
from MassTodonPy.Misc.cvxopt_wrapper import cvxopt_wrapper
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch, diag, incidence_matrix
from MassTodonPy.MatchMaker.CzMatch import CzMatch
from MassTodonPy.MatchMaker.RegularizedCzMatch import RegularizedCzMatch


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
# simple_matches.get_intensities()
# simple_matches.get_probabilities()
# simple_matches.get_branching_ratios()

matches = CzMatch(results, mol.precursor)
matches.get_intensities()
matches.get_probabilities()
matches.get_branching_ratios()

Q = matches.precursor.q
G = matches.graph.copy()


# lavish: all fragments lose cofragments
lavish = sum((Q - 1 - N.q) * I for N, I in G.nodes.data('intensity'))
node_intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
costs = matrix([float(ETnoD_PTR) for N,M, ETnoD_PTR in G.edges.data('ETnoD_PTR')])
edges_cnt = G.size()
equalities = incidence_matrix(G, len(node_intensities), edges_cnt)
inequalities = diag(-1.0, edges_cnt)
upper_bounds = matrix([0.0] * edges_cnt)

# print(costs, equalities, node_intensities, inequalities, upper_bounds)
primalstart = {}
primalstart['x'] = matrix([0.0] * edges_cnt)
primalstart['s'] = matrix([eps] * len(upper_bounds))

solution = solvers.conelp(c=costs,
                          G=inequalities,
                          h=upper_bounds,
                          A=equalities,
                          b=node_intensities,
                          primalstart=primalstart)

# Studying the max and min PTR and ETnoD on fragments

# removing self-loops
G.remove_edges_from([(N, M) for N, M in G.edges if N is M])
# removing trivial cases
G.remove_nodes_from([N for N, deg in G.degree if not deg])

if len(G) > 0:
    PTR_cnts = matrix([float(PTR) for N,M, PTR in G.edges.data('PTR')])
    ETnoD_cnts = matrix([float(ETnoD) for N,M, ETnoD in G.edges.data('ETnoD')])

    no_ETnoD = all(ETnoD == 0 for ETnoD in ETnoD_cnts)
    no_PTR = all(PTR == 0 for PTR in PTR_cnts)
    min_reaction_intensity = sum(I for N, M, I in G.edges.data('flow'))

    if no_ETnoD:
        PTR_min = PTR_max = min_reaction_intensity
        ETnoD_min = ETnoD_max = 0.0
    elif no_PTR:
        ETnoD_min = ETnoD_max = min_reaction_intensity
        PTR_min = PTR_max = 0.0
    else:
        node_intensities = matrix([float(I) for N, I in G.nodes.data('intensity')])
        edges_cnt = G.size()
        inequalities = sparse([incidence_matrix(G, len(node_intensities), edges_cnt),
                               diag(-1.0, edges_cnt)])
        inequality_constraints = matrix([node_intensities, matrix([0.0] * edges_cnt)])
        min_reaction_coef = matrix([float(ETnoD_PTR) for N,M, ETnoD_PTR in
                                    G.edges.data('ETnoD_PTR')]).T
        primalstart['x'] = matrix([0.0] * edges_cnt)
        primalstart['s'] = matrix([eps] * len(inequality_constraints))
        min_PTR = solvers.conelp(c=PTR_cnts,
                                 G=inequalities,
                                 h=inequality_constraints,
                                 A=min_reaction_coef,
                                 b=matrix(min_reaction_intensity),
                                 primalstart=primalstart)
        PTR_min = min_PTR['x']
        ETnoD_max = min_reaction_intensity - PTR_min
        max_PTR = solvers.conelp(c=-PTR_cnts,
                                 G=inequalities,
                                 h=inequality_constraints,
                                 A=min_reaction_coef,
                                 b=min_reaction_intensity,
                                 primalstart=primalstart)
        PTR_max = max_PTR['x']
        ETnoD_min = min_reaction_intensity - PTR_max

print(PTR_min, PTR_max, ETnoD_min, ETnoD_max)
