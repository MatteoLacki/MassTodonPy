from collections import defaultdict

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon

substanceP = get_dataset('substanceP')
modifications = defaultdict(dict)
for (number, group), mods in substanceP.precursor.modifications.items():
    modifications[number][group] = dict(mods)

precursor = {'name': 'substanceP',
             'fasta': substanceP.precursor.fasta,
             'charge': 3,
             'modifications': modifications}

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      precursor=precursor,
                      mz_tol=.05,
                      _devel=True)

masstodon.report.plot('/Users/matteo/Desktop/assigned_spectrum.html',
                      width= 1000)
