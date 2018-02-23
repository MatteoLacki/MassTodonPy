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

# agg_mols = {name: intensity for name, _, intensity in masstodon.report.aggregated_mols()}
# agg_mols
self = masstodon.report
# masstodon.report.M
# list(self.M.precursor.molecules())
# self.M.precursor._get_amino_acid(4, 'N')
self.M.cz_match.intensities['ETDorHTR_bond']
self.M.molecules




flatten_modification(self.M.precursor.modifications)


from collections import defaultdict, namedtuple

data_chunks = defaultdict(dict)
for mol in self.M.molecules:
    if mol.name is not 'precursor':
        mT, pos, cS = mol._molType_position_cleavageSite()
        data_chunks[cS][]


mol._molType_position_cleavageSite()



def iter_fragment_tags(precursor):
    for i, aa in enumerate(precursor.fasta):
        for precise_cleavage_site in precursor.groups:
            yield i, aa, precise_cleavage_site

list(iter_fragment_tags(self.M.precursor))


masstodon.report.aggregeted_fragment_intensities()

masstodon.write(path)


 from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_precursors import bokeh_aggregated_precursors
from MassTodonPy.Plot.bokeh_fragments_intensity import bokeh_fragments_intensity

bokeh_spectrum(masstodon, path + 'assigned_spectrum.html')
bokeh_aggregated_precursors(masstodon, path + 'aggregated_precusors.html')
bokeh_fragments_intensity(masstodon, path + 'fragment_intensities.html')

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet, FactorRange
from bokeh.transform import factor_cmap
from MassTodonPy.Misc.io import create_folder_if_needed

path3 = path + 'aggregeted_fragment_intensities.html'

create_folder_if_needed(path3)
output_file(path3)


afi = masstodon.report.aggregeted_fragment_intensities()


list(zip(afi['c_name'], afi['z_name']))


for c_name, z_name in zip(afi['c_name'], afi['z_name']):
    yield ('', c_name)


x = [ ('', ) for for _ in range(len(afi['c'])) ]

p = figure(x_range=)


p.vbar('x')
# Old plot
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
