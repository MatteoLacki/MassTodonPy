%load_ext autoreload
%autoreload 2

from collections import defaultdict
import csv

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Parsers.Paths import parse_path

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

masstodon.cz_match.probabilities
masstodon.cz_match.intensities

# Check the tests.

path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'

b = masstodon.report._bricks[0]
b.molecule.name

masstodon.report.plot(path=path1)
masstodon.report.write(path2)

_bricks = masstodon.report._bricks
_peak_groups = masstodon.report._peak_groups


_peak_groups
