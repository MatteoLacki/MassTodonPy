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

path="/Users/matteo/Desktop/assigned_spectrum.html"
plot = masstodon.report.plot(path=path,
                             height=1000,
                             width=2000)

plot_json = plot.to_json(include_defaults=True)

plot_json['renderers']

from bokeh.embed import components, autoload_static
from bokeh.resources import CDN

help(CDN)
help(autoload_static)
components(plot)
js, tag = autoload_static(plot, CDN, "/Users/matteo/Desktop/Bokeh")
