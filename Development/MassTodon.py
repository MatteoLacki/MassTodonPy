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
# masstodon.write(path)
self = masstodon.report
# from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
# from MassTodonPy.Plot.bokeh_aggregated_precursors import bokeh_aggregated_precursors
# from MassTodonPy.Plot.bokeh_fragments_intensity import bokeh_fragments_intensity
#
# bokeh_spectrum(masstodon, path + 'assigned_spectrum.html')
# bokeh_aggregated_precursors(masstodon, path + 'aggregated_precusors.html')
# bokeh_fragments_intensity(masstodon, path + 'fragment_intensities.html')
#


from MassTodonPy.Data.amino_acids import aa2name
from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet, FactorRange
from bokeh.transform import factor_cmap
from MassTodonPy.Misc.io import create_folder_if_needed


from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, Range1d, LabelSet, Label

offset = .2
data = defaultdict(list)
for d in self.iter_aggregated_fragments_intensity():
    for k, v in d.items():
        data[k].append(v)

data['x'] = list(range(len(masstodon.precursor.fasta)))
data['x_left'] = [i-offset for i in range(len(masstodon.precursor.fasta))]
data['x_right'] = [i+offset for i in range(len(masstodon.precursor.fasta))]


# aa = dict(x=[0, 0.5, 1, 1.5, 2],
#           y=[0, 0, 0, 0, 0],
#           names=['N', '-', 'C', '-', 'C'])

source = ColumnDataSource(data=data)
path3 = path + 'aggregeted_fragment_intensities.html'
create_folder_if_needed(path3)
output_file(path3)
p = figure()
p.xaxis.visible = False
left_bars = p.vbar(x='x_left', width=2*offset, top='left', color="orange", source=source)
right_bars= p.vbar(x='x_right', width=2*offset, top='right', color="navy", source=source)
left_labels = LabelSet(x='x_left', y='left', text='left_name', level='glyph', x_offset=-10,
                       source=source, render_mode='css')
right_labels = LabelSet(x='x_right', y='right', text='right_name', level='glyph', x_offset=-5,
                        source=source, render_mode='css')
aa_labels = LabelSet(x='x', y=0.0, text='aa', level='glyph', y_offset=-20, x_offset=-8,
                     source=source, render_mode='css')
p.add_layout(left_labels)
p.add_layout(right_labels)
p.add_layout(aa_labels)
hover_bars_l = HoverTool(renderers=[left_bars],
                         mode='vline',
                         tooltips=[('name', '@left_name'),
                                   ('intensity', "@left{0.}")])
p.add_tools(hover_bars_l)
hover_bars_r = HoverTool(renderers=[right_bars],
                         mode='vline',
                         tooltips=[('fragment', '@right_name'),
                                   ('intensity', "@right{0.}")])
p.add_tools(hover_bars_r)
show(p)


# second try
# from MassTodonPy.Misc.Iterators import alternate
# from bokeh.models import ColumnDataSource, FactorRange
#
# aas = data['aa']
# frags = ('c', 'z')
# data2 = dict(aa = data['aa'],
#              c  = data['left'],
#              z  = data['right'])
#
# def alternate_cz(C, Z, AA):
#     for c, z, aa in zip(C, Z, AA):
#         yield (c, aa)
#         yield (z, aa)
#
# x = list(alternate_cz(data['left_name'], data['right_name'], data['aa']))
# intensities = list(alternate(data['left'], data['right']))
#
# source2 = ColumnDataSource(data=dict(x=x, intensities=intensities))
# p = figure(x_range=FactorRange(*x))
# p.y_range.start = 0
# bars = p.vbar(x='x', width=2*offset, top='intensities', color="orange", source=source2, factors=aas)
#
# left_labels = LabelSet(x='x_left', y='left', text='left_name', level='glyph', x_offset=-15,
#                        source=source, render_mode='css')
# right_labels = LabelSet(x='x_right', y='right', text='right_name', level='glyph',
#                        source=source, render_mode='css')
# p.add_layout(left_labels)
# p.add_layout(right_labels)
# hover_bars_l = HoverTool(renderers=[left_bars],
#                          mode='vline',
#                          tooltips=[('name', '@left_name'),
#                                    ('intensity', "@left{0.}")])
# p.add_tools(hover_bars_l)
# hover_bars_r = HoverTool(renderers=[right_bars],
#                          mode='vline',
#                          tooltips=[('name', '@right_name'),
#                                    ('intensity', "@right{0.}")])
# p.add_tools(hover_bars_r)
# show(p)
