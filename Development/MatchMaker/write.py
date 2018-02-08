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

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications)

# Use this writer more often: only provide a flow of outcomes!!!


path = '/Users/matteo/Desktop/test/'
masstodon.write(path)
