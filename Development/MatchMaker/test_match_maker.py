from __future__ import absolute_import, division, print_function
from collections import Counter
import unittest

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch


def node_to_tuple(node):
    return node.type + str(node.no), node.q

mol = get_dataset('substanceP')
mols = list(mol.precursor.molecules())
mols_dict = {(m.name, m.q, m.g): m for m in mols}
results = [
    {'alphas': [
        {'estimate': 1000,
         'molecule': mols_dict[('c4', 1, 0)]},
        {'estimate': 5000,
         'molecule': mols_dict[('c5', 1, 0)]},
        {'estimate': 8000,
         'molecule': mols_dict[('c6', 1, 0)]}],
     'status': 'optimal'},
    {'alphas': [
        {'estimate': 1200,
         'molecule': mols_dict[('z7', 1, 0)]},
        {'estimate': 4800,
         'molecule': mols_dict[('z6', 1, 0)]},
        {'estimate': 8300,
         'molecule': mols_dict[('z5', 1, 0)]}],
     'status': 'optimal'},
    {'alphas': [
        {'estimate': 40000,
         'molecule': mols_dict[('precursor', 3, 0)]},
        {'estimate': 30000,
         'molecule': mols_dict[('precursor', 2, 0)]},
        {'estimate': 20000,
         'molecule': mols_dict[('precursor', 1, 0)]}],
     'status': 'optimal'}]

matches = SimpleCzMatch(results, mol.precursor)

expected_nodes = set([('c4', 1), ('z7', 1), ('c5', 1),
                      ('z6', 1), ('c6', 1), ('z5', 1)])
obtained_nodes = set(node_to_tuple(n)
                     for n in matches.graph.nodes)

expected_nodes
obtained_nodes

expected_edges = set(frozenset(e) for e in
                     [[('c4', 1), ('z7', 1)],
                      [('c5', 1), ('z6', 1)],
                      [('c6', 1), ('z5', 1)],
                      [('c4', 1), ('c4', 1)],
                      [('c5', 1), ('c5', 1)],
                      [('c6', 1), ('c6', 1)],
                      [('z5', 1), ('z5', 1)],
                      [('z6', 1), ('z6', 1)],
                      [('z7', 1), ('z7', 1)]])



obtained_edges = set(frozenset([node_to_tuple(n),
                                node_to_tuple(m)])
                     for n, m in matches.graph.edges)
expected_edges
obtained_edges
