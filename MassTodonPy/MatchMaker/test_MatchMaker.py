from __future__ import absolute_import, division, print_function
from collections import Counter
import unittest

from MassTodonPy.Data.get_dataset import get_dataset
# from MassTodonPy.MatchMaker.match_cz_ions import czMatchMakerBasic
from MassTodonPy.MatchMaker.Base_cz_match import Base_cz_match
from MassTodonPy.MatchMaker.match_cz_ions import czMatchMakerIntermediate
from MassTodonPy.MatchMaker.match_cz_ions import czMatchMakerAdvanced

# TODO add test for the optimality of the fragment matching

def node_to_tuple(node):
    return node.type + str(node.no), node.q

class TestPeakPicker(unittest.TestCase):
    def setUp(self):
        """Set up a method."""

        self.mol = get_dataset('substanceP')
        mols = list(self.mol.precursor.molecules())
        mols_dict = {(m.name, m.q, m.g): m for m in mols}
        self.results = [
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


    def tearDown(self):
        """Tear down a method."""
        pass

    def test_basic_matchmaker(self):
        print("Testing basic matchmaker: graph creation.")

        matches = Base_cz_match(self.results, self.mol.precursor)
        expected_nodes = set([('c4', 1), ('z7', 1), ('c5', 1),
                              ('z6', 1), ('c6', 1), ('z5', 1)])
        obtained_nodes = set(node_to_tuple(n) for n in matches.graph.nodes)
        self.assertEqual(expected_nodes, obtained_nodes)

        expected_edges = set([frozenset([('c4', 1), ('z7', 1)]),
                              frozenset([('c5', 1), ('z6', 1)]),
                              frozenset([('c6', 1), ('z5', 1)])])
        obtained_edges = set(frozenset([node_to_tuple(n), node_to_tuple(m)])
                             for n, m in matches.graph.edges)
        
        self.assertEqual(expected_edges, obtained_edges)

    def test_intermediate_matchmaker(self):
        print("Testing basic matchmaker: graph creation.")

        matches = czMatchMakerIntermediate(self.results, self.mol.precursor)
        expected_nodes = set([('c4', 1, 0), ('z7', 1, 0), ('c5', 1, 0),
                              ('z6', 1, 0), ('c6', 1, 0), ('z5', 1, 0)])
        obtained_nodes = set(matches.graph.nodes)
        self.assertEqual(expected_nodes, obtained_nodes)
        expected_edges = set([frozenset([('c4', 1, 0), ('z7', 1, 0)]),
                              frozenset([('c5', 1, 0), ('z6', 1, 0)]),
                              frozenset([('c6', 1, 0), ('z5', 1, 0)])])
        obtained_edges = set(frozenset(e) for e in matches.graph.edges)
        self.assertEqual(expected_edges, obtained_edges)

    def test_advanced_matchmaker(self):
        print("Testing advanced matchmaker: graph creation.")
        pass

if __name__ == "__main__":
    unittest.main()
