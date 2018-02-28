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


masstodon = MassTodon(spectrum=substanceP.spectrum,
                      min_intensity=10.0,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications,
                      _verbose=True)

from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_fragments import bokeh_aggregated_fragments
from MassTodonPy.Plot.bokeh_aggregated_fragments_estimated import bokeh_aggregated_fragments_estimated

bokeh_spectrum(masstodon, show=True)
bokeh_aggregated_fragments(masstodon, show=True)
bokeh_aggregated_fragments_estimated(masstodon, show=True)
