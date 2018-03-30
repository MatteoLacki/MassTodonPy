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

print(modifications)

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      min_intensity=10.0,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications,
                      _verbose=False)

from MassTodonPy.Plot import bokeh_spectrum
from MassTodonPy.Plot import bokeh_aggregated_precursors
from MassTodonPy.Plot import bokeh_aggregated_fragments
from MassTodonPy.Plot import bokeh_aggregated_fragments_estimated

bokeh_spectrum(masstodon, path=path + 'assignes_spectrum.html', show=True)
bokeh_aggregated_fragments(masstodon, show=True)
bokeh_aggregated_fragments_estimated(masstodon, show=True)
bokeh_aggregated_precursors(masstodon, show=True)
