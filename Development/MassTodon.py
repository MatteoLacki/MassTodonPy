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
masstodon.write(path)

from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_precursors import bokeh_aggregated_precursors
from MassTodonPy.Plot.bokeh_fragments_intensity import bokeh_fragments_intensity

bokeh_spectrum(masstodon, path + 'assigned_spectrum.html', width=300, height=400)
bokeh_aggregated_precursors(masstodon, path + 'aggregated_precusors.html')
bokeh_fragments_intensity(masstodon, path + 'fragment_intensities.html')



from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet
from MassTodonPy.Misc.os import create_folder_if_needed

path3 = path + 'aggregeted_fragment_intensities.html'


create_folder_if_needed(path3)
output_file(path3)

afi = masstodon.report.aggregeted_fragment_intensities()
afi['z_minus'] = [-v for v in afi['z']]

# p = figure(plot_width=1000,
#            plot_height=400)
p = figure()
p.xaxis.axis_label = 'Cleavage Site'
p.yaxis.axis_label = 'Estimated Intensity'
afi['x'] = list(range(len(afi['z'])))
data = ColumnDataSource(afi)
bars_c = p.vbar(source=data, x='x', top='c', width=.8)
hover_bars_c = HoverTool(renderers=[bars_c],
                         tooltips=[('name', '@c_name'),
                                   ('intensity', "@c{0,0}")])
p.add_tools(hover_bars_c)
bars_z = p.vbar(source=data, x='x', top='z_minus', width=.8, color='red')
hover_bars_z = HoverTool(renderers=[bars_z],
                         tooltips=[('name', '@z_name'),
                                   ('intensity', "@z{0,0}")])
p.add_tools(hover_bars_z)
# labels_c = LabelSet(source=data, x='x', y='c')
# p.add_layout(labels_c)
# labels_z = LabelSet(source=data, x='x', y='z')
# p.add_layout([labels_z, labels_c])
show(p)
