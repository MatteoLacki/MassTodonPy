"""Testing test_formulator."""

import unittest
from MassTodonPy.Formulator.formulator import make_formulas
from MassTodonPy import get_data

class TestFormulator(unittest.TestCase):
    def setUp(self):
        """Set up a method."""
        pass

    def tearDown(self):
        """Tear down a method."""
        pass

    def test_make_formulas(self):
        print("Testing make_formulas")
        expected_molecules = [('precursor', 'C63H98N18O13S1', 11, 1, 0), ('precursor', 'C63H98N18O13S1', 11, 1, 1), ('precursor', 'C63H98N18O13S1', 11, 2, 0), ('c2', 'C11H22N6O2', 2, 1, -1), ('c2', 'C11H22N6O2', 2, 1, 0), ('c4', 'C22H41N9O4', 4, 1, -1), ('c4', 'C22H41N9O4', 4, 1, 0), ('c5', 'C27H49N11O6', 5, 1, -1), ('c5', 'C27H49N11O6', 5, 1, 0), ('c6', 'C32H57N13O8', 6, 1, -1), ('c6', 'C32H57N13O8', 6, 1, 0), ('c7', 'C41H66N14O9', 7, 1, -1), ('c7', 'C41H66N14O9', 7, 1, 0), ('c8', 'C50H75N15O10', 8, 1, -1), ('c8', 'C50H75N15O10', 8, 1, 0), ('c9', 'C52H78N16O11', 9, 1, -1), ('c9', 'C52H78N16O11', 9, 1, 0), ('c10', 'C58H89N17O12', 10, 1, -1), ('c10', 'C58H89N17O12', 10, 1, 0), ('z1', 'C5H10N1O1S1', 1, 1, 0), ('z1', 'C5H10N1O1S1', 1, 1, 1), ('z2', 'C11H21N2O2S1', 2, 1, 0), ('z2', 'C11H21N2O2S1', 2, 1, 1), ('z3', 'C13H24N3O3S1', 3, 1, 0), ('z3', 'C13H24N3O3S1', 3, 1, 1), ('z4', 'C22H33N4O4S1', 4, 1, 0), ('z4', 'C22H33N4O4S1', 4, 1, 1), ('z5', 'C31H42N5O5S1', 5, 1, 0), ('z5', 'C31H42N5O5S1', 5, 1, 1), ('z6', 'C36H50N7O7S1', 6, 1, 0), ('z6', 'C36H50N7O7S1', 6, 1, 1), ('z7', 'C41H58N9O9S1', 7, 1, 0), ('z7', 'C41H58N9O9S1', 7, 1, 1), ('z9', 'C52H77N12O11S1', 9, 1, 0), ('z9', 'C52H77N12O11S1', 9, 1, 1), ('z11', 'C63H96N17O13S1', 11, 1, 0), ('z11', 'C63H96N17O13S1', 11, 1, 1)]

        mol = get_data('substanceP')
        frags = make_formulas(fasta=mol["fasta"],
                              Q=2,
                              modifications=mol["modifications"])

        molecules = list(frags.makeMolecules())
        self.assertEqual(molecules, expected_molecules)


if __name__ == "__main__":
    unittest.main()
