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

# Use this writer more often: only provide a flow of outcomes!!!

def write_from_buffer(things, path, name):
    file_path, file_name, file_ext = parse_path(path)
    assert file_ext in ('.csv', '.tsv'), "Writing only to csv or tsv."
    delimiter = ',' if file_ext == '.csv' else '\t'
    path_things = path.replace('.', '_intensities.')
    with open(path_things, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=delimiter)
        for res in things:
            writer.writerow(res)

path = '/Users/matteo/Desktop/pairing.csv'
write_from_buffer(masstodon.simple_cz_match._iter_intensities(), path, '_intensities')

masstodon.cz_match.intensities
masstodon.cz_match.probabilities


path_probs = '/Users/matteo/Desktop/pairing.csv'
write_from_buffer(masstodon.simple_cz_match._iter_probabilities(), path_probs, '_probs')
