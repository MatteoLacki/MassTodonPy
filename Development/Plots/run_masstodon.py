from collections import Counter, defaultdict

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Plotting.plot_buffers import buffers
from MassTodonPy.Misc.sorting import sort_by_first
from MassTodonPy.Spectra.Spectrum import Spectrum

def get_masstodon_results():
    substanceP = get_dataset('substanceP')
    # substanceP.spectrum.plot()

    modifications = defaultdict(dict)
    for (number, group), mods in substanceP.precursor.modifications.items():
        modifications[number][group] = dict(mods)

    precursor = {'name': 'substanceP',
                 'fasta': substanceP.precursor.fasta,
                 'charge': 3,
                 'modifications': modifications}

    masstodon = MassTodon(spectrum=substanceP.spectrum,
                          precursor=precursor,
                          mz_precision=.05,
                          _devel=True)
    return masstodon
