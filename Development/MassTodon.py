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
                      min_intensity=10.0,
                      mz_tol=.05,
                      fasta=substanceP.precursor.fasta,
                      charge=3,
                      name='substanceP',
                      modifications=modifications)

path = '/Users/matteo/Desktop/test/'
path2= path + 'assigned_spectrum.csv'

from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_precursors import plot_aggregated_precursors

bokeh_spectrum(masstodon, path + 'assigned_spectrum.html')
plot_aggregated_precursors(masstodon, path + 'aggregated_precusors.html')








from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet
from MassTodonPy.Misc.os import create_folder_if_needed


create_folder_if_needed(path3)
output_file(path3)

intensities = masstodon.report.get_aggregated_precursors()

charges = list('q = {0}'.format(x) for x in range(1, 1+len(intensities)))
p = figure(plot_width=400,
           plot_height=400,
           x_range=charges)
p.xaxis.axis_label = 'Charge'
p.yaxis.axis_label ='Estimated Intensity'

data = ColumnDataSource({'intensity': intensities,
                         'charge': charges})
bars = p.vbar(source=data, x='charge', top='intensity', width=1)
hover_bars = HoverTool(renderers=[bars],
                       tooltips=[('charge', '@charge'),
                                 ('intensity', "@intensity{0,0}")],
                       mode='vline')
p.add_tools(hover_bars)
show(p)
