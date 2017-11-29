%load_ext autoreload
%autoreload 2

# simply take a few molecules with different quenched charge
# and make the spectrum equal to their isotopic distributions.

# which is why it would be nice to have a scalar multiplication
# of the measure first.

# The test should check what m/z are present and where they should be

from collections import Counter
from networkx import connected_component_subgraphs as components
import matplotlib.pyplot as plt
import networkx as nx

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.PeakPicker.PeakPicker2 import get_deconvolution_problems
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum

options = {'node_color': 'black',
           'node_size': 2,
           'width': 1,
           'with_labels': True,
           'font_size': 10}


subP = get_dataset('substanceP')

precursors = list(mol for mol in subP.precursor.molecules()
                  if mol.name is 'precursor')

# for mol in precursors:
#     spec = mol.isotopologues()
#     spec.round_mz(2)
    # print(spec.mz)

spectrum = sum(mol.isotopologues() for mol in precursors)
spectrum = ExperimentalSpectrum(mz=spectrum.mz,
                                intensity=100000 * spectrum.probability)
spectrum.round_mz(precision=2)

DG = get_deconvolution_problems(precursors,
                                spectrum,
                                mz_tol=.05,
                                mz_precision=2)

DGs = list(DG)

stats = [Counter(N[0] for N in DG) for DG in DGs ]
stats = set([(s['M'], s['I'], s['G'])for s in stats])



nx.draw(DGs[0], **options)
plt.figure(figsize=(12,12))
plt.show()

nx.draw(DGs[1], **options)
plt.figure(figsize=(12,12))
plt.show()

nx.draw(DGs[2], **options)
plt.figure(figsize=(12,12))
plt.show()

# Counter(N[0] for N in DG)



# for N in DG:
#     if N[0] is 'G':
#         print(DG[N])
#
#
#
# # Maybe we poorly number G nodes?
# Counter((a[0], b[0]) for a, b in DG.edges)
#
#
# DGs[0].nodes(data=True)
#
#
#
# DGs = list(components(DG))
# len(DGs)
