"""Testing the setting up of the deconvolution problems."""
from __future__ import absolute_import, division, print_function
from collections import Counter
import unittest

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Deconvolution.Deconvolve import deconvolve
from MassTodonPy.Spectra.Spectrum import Spectrum


class TestPeakPicking(unittest.TestCase):
    def test_deconvolve(self):
        print("Testing the get_deconvolution_problems function.")
        # R_ = real
        R_stats = {(1, 11, 7), (2, 27, 8), (3, 54, 9)}
        subP = get_dataset('substanceP')
        precursors = list(m for m in subP.precursor.molecules()
                          if m.name is 'precursor')
        spectrum = sum(m.isotopologues() for m in precursors)
        spectrum = Spectrum(mz=spectrum.mz,
                            intensity=100000 * spectrum.probability)
        spectrum.round_mz(precision=2)
        DGs = deconvolve(precursors,
                         spectrum,
                         'Matteo',
                         mz_tol=.05,
                         _isospec_args={'mz_precision': 2})
        # E_ = expected
        E_stats = [Counter(N[0] for N in DG) for DG in DGs ]
        E_stats = set([(s['M'], s['I'], s['G']) for s in E_stats])
        self.assertEqual(R_stats, E_stats)

if __name__ == "__main__":
    unittest.main()
