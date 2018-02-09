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
                      min_intensity=100.0,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications)

path = '/Users/matteo/Desktop/test/'
path1= path + 'assigned_spectrum.html'
path2= path + 'assigned_spectrum.csv'
# plot = masstodon.plot(path1, width= 1000)
masstodon.write(path)
# masstodon.report._spectrum.total_intensity()
# masstodon.report._spectrum.low_spectrum.total_intensity()
# masstodon.report._spectrum.low_spectrum.l1()
# masstodon.report._spectrum.low_spectrum.l2()
# list(masstodon.report.iter_global_quality_fits())
# masstodon.report.global_quality_fits_stats()

for sol in masstodon._solutions:
    print(sol.global_fit_quality())

masstodon.report._spectrum.low_spectrum.l1

%load_ext autoreload
%autoreload 2

from MassTodonPy.Spectra.Spectrum import Spectrum

s = Spectrum()
s.total_intensity()
s.l1()
s.l2()
