from __future__ import absolute_import, division, print_function
from collections import defaultdict, Counter
import unittest

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon


class TestMassTodonOnSubstanceP(unittest.TestCase):
    def setUp(self):
        """Load data on Substance P."""

        substanceP = get_dataset('substanceP')
        modifications = defaultdict(dict)
        for (number, group), mods in substanceP.precursor.modifications.items():
            modifications[number][group] = dict(mods)
        
        self.args = {'spectrum':        substanceP.spectrum,
                     'min_intensity':   50.0,
                     'fasta':           substanceP.precursor.fasta,
                     'name':            'substanceP',
                     'modifications':   modifications,
                     'charge':          3,
                     'mz_tol':          .05}
        
    def run_masstodon(self):
        """Run one session of the MassTodon on the Substance P data."""
        masstodon = MassTodon(**self.args)

    def test_some_times(self):
        """Test the MassTodon on substance P, M times."""
        M = 100
        res = Counter() 
        for i in range(M):
            try:
                masstodon = self.run_masstodon()
                res['success'] += 1
            except ValueError as e:
                res['failure'] += 1
        
        self.assertEqual(res['success'], M)
        self.assertEqual(res['failure'], 0)
                                                  
#class TestPeakPicker(unittest.TestCase):
#    def setUp(self):
#        """Set up a method."""
#
#        mol = get_dataset('substanceP')
#        mols = list(mol.precursor.molecules())
#        mols_dict = {(m.name, m.q, m.g): m for m in mols}
#        results = [
#            {'alphas': [
#                {'estimate': 1000,
#                 'molecule': mols_dict[('c4', 1, 0)]},
#                {'estimate': 5000,
#                 'molecule': mols_dict[('c5', 1, 0)]},
#                {'estimate': 8000,
#                 'molecule': mols_dict[('c6', 1, 0)]}],
#             'status': 'optimal'},
#            {'alphas': [
#                {'estimate': 1200,
#                 'molecule': mols_dict[('z7', 1, 0)]},
#                {'estimate': 4800,
#                 'molecule': mols_dict[('z6', 1, 0)]},
#                {'estimate': 8300,
#                 'molecule': mols_dict[('z5', 1, 0)]}],
#             'status': 'optimal'},
#            {'alphas': [
#                {'estimate': 40000,
#                 'molecule': mols_dict[('precursor', 3, 0)]},
#                {'estimate': 30000,
#                 'molecule': mols_dict[('precursor', 2, 0)]},
#                {'estimate': 20000,
#                 'molecule': mols_dict[('precursor', 1, 0)]}],
#             'status': 'optimal'}]
#
#        molecules = [m['molecule'] for res in results for m in res['alphas']]
#        # for res in results:
#        #     for x in res['alphas']:
#        #         mol = VoidClass()
#        #         mol.intensity = x['estimate']
#        #         x['molecule'].intensity =
#        #         brick = Brick(molecule=x['molecule'])
#        #         bricks.append(brick)
#        #
#        # # Preparing a phoney MassTodon
#        # self.precursor = VoidClass()
#        self.precursor_charge = 3
#        self.molecules = molecules
#
#    def tearDown(self):
#        """Tear down a method."""
#        pass
#
#    def test_SimpleCzMatch(self):
#        print("Testing SimpleCzMatch graph creation.")
#
#        # precursor_charge = 3
#        matches = SimpleCzMatch(molecules=self.molecules,
#                                precursor_charge=self.precursor_charge)
#
#        expected_nodes = set([('c4', 1), ('z7', 1), ('c5', 1),
#                              ('z6', 1), ('c6', 1), ('z5', 1)])
#        obtained_nodes = set(node_to_tuple(n)
#                             for n in matches.graph.nodes)
#        self.assertEqual(expected_nodes, obtained_nodes)
#
#        expected_edges = set(frozenset(e) for e in
#                     [[('c4', 1), ('z7', 1)],
#                      [('c5', 1), ('z6', 1)],
#                      [('c6', 1), ('z5', 1)],
#                      [('c4', 1), ('c4', 1)],
#                      [('c5', 1), ('c5', 1)],
#                      [('c6', 1), ('c6', 1)],
#                      [('z5', 1), ('z5', 1)],
#                      [('z6', 1), ('z6', 1)],
#                      [('z7', 1), ('z7', 1)]])
#
#        obtained_edges = set(frozenset([node_to_tuple(n),
#                                        node_to_tuple(m)])
#                             for n, m in matches.graph.edges)
#        self.assertEqual(expected_edges, obtained_edges)
#
    # def test_intermediate_matchmaker(self):
    #     print("Testing basic matchmaker: graph creation.")
    #     matches = czMatchMakerIntermediate(masstodon=self.masstodon)

    #     expected_nodes = set([('c4', 1, 0), ('z7', 1, 0), ('c5', 1, 0),
    #                           ('z6', 1, 0), ('c6', 1, 0), ('z5', 1, 0)])
    #     obtained_nodes = set(matches.graph.nodes)
    #     self.assertEqual(expected_nodes, obtained_nodes)

    #     expected_edges = set(frozenset(e) for e in
    #                         [[('c4', 1, 0), ('z7', 1, 0)],
    #                          [('c5', 1, 0), ('z6', 1, 0)],
    #                          [('c6', 1, 0), ('z5', 1, 0)]])
    #     obtained_edges = set(frozenset(e) for e in matches.graph.edges)
    #     self.assertEqual(expected_edges, obtained_edges)

    # def test_advanced_matchmaker(self):
    #     print("Testing advanced matchmaker: graph creation.")
    #     pass

if __name__ == "__main__":
    unittest.main()
