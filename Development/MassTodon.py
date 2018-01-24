%load_ext autoreload
%autoreload 2

from collections import defaultdict

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon

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

masstodon.reporter.plot(file_path='/Users/matteo/Desktop/assigned_spectrum.html')
