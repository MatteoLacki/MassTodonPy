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
precursor = {'name': 'substanceP',
             'fasta': substanceP.precursor.fasta,
             'charge': 3,
             'modifications': modifications}
masstodon = MassTodon(spectrum=substanceP.spectrum,
                      precursor=precursor,
                      mz_precision=.05,
                      simple_cz_match=True,
                      _devel=True)


masstodon.simple_cz_match.probabilities
masstodon.cz_match.probabilities


simple = masstodon.simple_cz_match.intensities
sum(v for k, v in simple['ETD_bond'].items())


adv = masstodon.cz_match.intensities
simple['ETD_bond']
adv['ETD_bond']

for k in simple:
    print("".format(k, simple[k], adv[k]))


from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch

masstodon.molecules


simple_cz_match = SimpleCzMatch(masstodon.molecules,
                                masstodon.precursor.q)

simple_cz_match.probabilities
simple_cz_match.intensities
# masstodon.cz_match.probabilities
# masstodon.cz_match.intensities

# Check the tests.

.path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'

masstodon.report.plot(path1)
masstodon.report.write(path2)

_bricks = masstodon.report._bricks
_peak_groups = masstodon.report._peak_groups
_clusters = masstodon.report._clusters
