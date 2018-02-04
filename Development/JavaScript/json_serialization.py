%load_ext autoreload
%autoreload 2


from collections import defaultdict
import csv

from bokeh.core.json_encoder import serialize_json, BokehJSONEncoder

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

path = '/Users/matteo/Desktop/'
path1= path + 'assigned_spectrum.html'

plot = masstodon.report.plot(path1, width= 1000)

help(BokehJSONEncoder)

from bokeh.document import Document

doc = Document()
doc.add_root(plot)


from bokeh.resources import CDN
from bokeh.embed import autoload_static

help(plot)
plot.to_json(include_defaults=True)

help(autoload_static)
js, tag = autoload_static(plot, CDN, "some/path")


help(Document)

serialize_json(plot)
