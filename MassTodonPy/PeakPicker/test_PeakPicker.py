"""Testing the setting up of the deconvolution problems."""
from __future__ import absolute_import, division, print_function
from collections import Counter
import unittest

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.PeakPicker.PeakPicker import get_deconvolution_problems
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum


class TestPeakPicker(unittest.TestCase):
    def setUp(self):
        """Set up a method."""
        pass

    def tearDown(self):
        """Tear down a method."""
        pass

    def test_get_deconvolution_problems(self):
        print("Testing the get_deconvolution_problems function.")

        real_stats = {(1, 11, 7), (2, 27, 8), (3, 54, 9)}

        subP = get_dataset('substanceP')
        precursors = list(mol for mol in subP.precursor.molecules()
                          if mol.name is 'precursor')

        spectrum = sum(mol.isotopologues() for mol in precursors)
        spectrum = ExperimentalSpectrum(mz=spectrum.mz,
                                        intensity=100000 * spectrum.probability)
        spectrum.round_mz(precision=2)

        DGs = get_deconvolution_problems(precursors,
                                         spectrum,
                                         mz_tol=.05,
                                         mz_precision=2)

        expected_stats = [Counter(N[0] for N in DG) for DG in DGs ]
        expected_stats = set([(s['M'], s['I'], s['G'])
                              for s in expected_stats])

        self.assertEqual(real_stats, expected_stats)


if __name__ == "__main__":
    unittest.main()
