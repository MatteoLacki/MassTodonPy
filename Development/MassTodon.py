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

# substanceP.spectrum.plot('plot.html', width=1000, height=500)

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications)

path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'
# plot = masstodon.report.plot(path1, width= 1000)
path = '/Users/matteo/Desktop/test/'
masstodon.write(path)


# masstodon.report.to_json(path + 'json_4_plot/test.json')
sol = masstodon._solutions[0]




masstodon.molecules

sol.report()
