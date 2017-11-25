%load_ext autoreload
%autoreload 2

from MassTodonPy.Spectra.Measure import Measure
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum
from MassTodonPy.Data.get_dataset import get_dataset

subP = get_dataset('substanceP')

measure_A = ExperimentalSpectrum(mz=subP.spectrum.mz,
                                 intensity=subP.spectrum.intensity)
measure_B = ExperimentalSpectrum(mz=subP.spectrum.mz,
                                 intensity=subP.spectrum.intensity)

for atom, mass in measure_A:
    print(atom, mass)

measure_C = measure_A + measure_B
measure_C

measure_A.intensity[0]
measure_B.intensity[0]
measure_C.intensity[0]

measure_A += measure_B
measure_A.intensity[0]

measure_D = sum([measure_A, measure_B, measure_C])
measure_D.intensity[0]

measure_E = ExperimentalSpectrum(mz=subP.spectrum.mz,
                                 intensity=subP.spectrum.intensity)

measure_D + measure_E


spec_A = ExperimentalSpectrum(mz=subP.spectrum.mz,
                              intensity=subP.spectrum.intensity)

spec_A
spec_A.round_atoms(2)

spec_A

spec_A.trim(20.0)
spec_A


spec_A.split_measure(100)
spec_A

spec_A.get_cut_off_value(.99)
spec_A

len(spec_A)


spec = ExperimentalSpectrum(*subP.spectrum)

spec.mz

list(spec[1800, 2900])
