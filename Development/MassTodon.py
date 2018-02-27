%load_ext autoreload
%autoreload 2


from collections import defaultdict
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Parsers.Paths import parse_path

substanceP = get_dataset('substanceP')
modifications = defaultdict(dict)
for (number, group), mods in substanceP.precursor.modifications.items():
    modifications[number][group] = dict(mods)
modifications = dict(modifications)
# substanceP.spectrum.plot('plot.html', width=1000, height=500)

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      min_intensity=10.0,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications,
                      _verbose=True)

# {(10, 'C_carbo'): Formula('H': 1, 'O': -1, 'N': 1)}
# {(9, 'C_carbo'): Formula('H': 1, 'O': -1, 'N': 1)}
# Counter({'I': 288, 'E': 61, 'M': 33})

masstodon.molecules
masstodon.spectrum.intensity.sum()

# test_args = {'_L1_flow': 0.01,
#  '_L1_intensity': 0.01,
#  '_L2_intensity': 0.01,
#  '_max_times': 10,
#  '_maxiters': 1000,
#  'block_prolines': True,
#  'blocked_fragments': set(['c0']),
#  'charge': 3,
#  'deconvolution_method': 'Matteo',
#  'distance_charges': 5,
#  'fasta': 'RPKPQQFFGLM',
#  'fragments': 'cz',
#  'joint_probability': 0.999,
#  'max_buffer_len': 0.5,
#  'min_intensity': 2.2204460492503131e-16,
#  'min_prob_per_molecule': 0.7,
#  'modifications': {10: {'C_carbo': {'H': 1, 'N': 1, 'O': -1}}},
#  'mz_tol': 0.05,
#  'name': 'SubstanceP',
#  'ni2': 0.1,
#  'output': 'output/',
#  'percent_top_peaks': 1.0,
#  'show_progress': False,
#  'sigma2': 0.1,
#  'spectrum': '/Users/matteo/Documents/MassTodon/Tests/CLI/subP_spectrum.txt',
#  'verbose': True}
#
# masstodon = MassTodon(**test_args)
path = '/Users/matteo/Desktop/test/'

from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
bokeh_spectrum(masstodon, show=True)

self = masstodon.report

from MassTodonPy.Plot.bokeh_aggregated_fragments import bokeh_aggregated_fragments
from MassTodonPy.Plot.bokeh_aggregated_fragments_estimated import bokeh_aggregated_fragments_estimated

bokeh_aggregated_fragments(masstodon, show=True)
bokeh_aggregated_fragments_estimated(masstodon, show=True)
