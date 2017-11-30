%load_ext autoreload
%autoreload 2

from __future__ import absolute_import, division, print_function
import numpy as np

from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator
from MassTodonPy.Misc.binomial import binomial

iso_calc = IsotopeCalculator(mz_precision=0,
                             _isotope_masses={'T': [1,1000]},
                             _isotope_probabilities={'T': [0.5, 0.5]})

distribution = iso_calc.get_envelope(formula='T10',
                                     joint_probability=1.0)

expected_atoms = np.array([1.*(10-i) + 1000.*i for i in range(11)])
expected_masses = np.array([binomial(10,i) for i in range(11)]) / 2**10

all(expected_atoms == distribution.atoms)
all(np.abs(expected_masses - distribution.masses) < 1e-12)
