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
                      mz_tol=.05,
                      _devel=True)

simple = masstodon.simple_cz_match.intensities
adv = masstodon.cz_match.intensities

path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'

masstodon.report.plot(path1, width= 1000)
masstodon.report.write(path2)
# Add an aggregated report based on the list of references of molecules.
# Order by estimated intensity.
# Can inplace reordering screw anything? Think about it!

masstodon.molecules.sort(key=lambda m: m.intensity, reverse=True)
m = masstodon.molecules[0]
m.name
m.source.fasta

m.source.formula.str_with_charges(3, 0)
m.intensity
