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

negative_g_mols = [mol for mol in masstodon.molecules if mol.intensity > 0.0 and mol.g < 0]


simple = masstodon.simple_cz_match.intensities
adv = masstodon.cz_match.intensities

simple['ETDorHTR_bond']
adv['ETDorHTR_bond']

# There must be vastly
for k in simple['ETDorHTR_bond']:
    print("{0}\t{1}".format(k,
                            adv['ETDorHTR_bond'][k] - simple['ETDorHTR_bond'][k]))

from MassTodonPy.MatchMaker.SimpleCzMatch import SimpleCzMatch

masstodon.molecules


simple_cz_match = SimpleCzMatch(masstodon.molecules,
                                masstodon.precursor.q)

simple_cz_match.probabilities
simple_cz_match.intensities
