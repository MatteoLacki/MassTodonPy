%load_ext autoreload
%autoreload 2


from collections import defaultdict
import csv

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

simple = masstodon.simple_cz_match.intensities
adv = masstodon.cz_match.intensities

path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'

masstodon.report.plot(path1, width= 1000)
masstodon.report.write(path2)
