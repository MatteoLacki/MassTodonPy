%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from collections import Counter, defaultdict
import numpy as np

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Misc.plot_buffers import buffers
from MassTodonPy.Misc.sorting import sort_by_first
from MassTodonPy.Spectra.Spectrum import Spectrum

substanceP = get_dataset('substanceP')
# substanceP.spectrum.plot()

modifications = defaultdict(dict)
for (number, group), mods in substanceP.precursor.modifications.items():
    modifications[number][group] = dict(mods)

precursor = {'name': 'substanceP',
             'fasta': substanceP.precursor.fasta,
             'charge': 3,
             'modifications': modifications}

masstodon = MassTodon(spectrum=substanceP.spectrum,
                      precursor=precursor,
                      mz_precision=.05,
                      _devel=True)

# sol = masstodon._solutions[0]
# sol
# sol.nodes['M1']
# sol.nodes['G1']
# sol.nodes['I1']
# sol['G1']['I1']
# sol['G1']['I19']

#TODO extract M information.
def data(solutions, mz_digits, verbosity=0):
    base_width = 10**(-mz_digits)
    for sol in solutions:
        for G in sol:
            if G[0] is 'G':
                try:
                    min_mz = sol.node[G]['min_mz'] - base_width/2
                    max_mz = sol.node[G]['max_mz'] + base_width/2
                except KeyError:
                    mz = sol.node[G]['mz']
                    min_mz = mz - base_width/2
                    max_mz = mz + base_width/2
                intensity = sol.node[G]['intensity']
                if intensity > 0 or estimate > 0:
                    intensity_d = intensity/(max_mz - min_mz)
                    estimate = sol.node[G]['estimate']
                    estimate_d = estimate/(max_mz - min_mz)
                    out_counter = Counter({'intensity': intensity,
                                           'intensity_d': intensity_d,
                                           'estimate': estimate,
                                           'estimate_d': estimate_d})
                    if verbosity > 0:
                        out_formulas = []
                        # What is going on now with M - I1 - G
                        #                           M - I2 - G ? They are glued.
                        # We mostly cannot tell apart different isotopologues.
                        for I in sol[G]:
                            GM_estimate = sol[G][I]['estimate']
                            for M in sol[I]:
                                if M[0] is 'M':
                                    mol = sol.node[M]['molecule']
                                    out_formulas.append((str(mol.formula),
                                                         mol.q,
                                                         mol.g,
                                                         GM_estimate))
                        yield (min_mz, max_mz), (out_counter, out_formulas)
                    else:
                        yield (min_mz, max_mz), out_counter


# D = defaultdict(Counter)
# for interval, data in data(masstodon._solutions, masstodon.mz_digits):
#     D[interval] += data

D = defaultdict(Counter)
E = defaultdict(list)
for interval, (data, formulas) in data(masstodon._solutions,
                                       masstodon.mz_digits,
                                       1):
    D[interval] += data
    E[interval].extend(formulas)

mz_L, mz_R, I, I_d, E, E_d = list(zip(*sorted((l, r,
                                               c['intensity'],
                                               c['intensity_d'],
                                               c['estimate'],
                                               c['estimate_d'])
                                              for (l, r), c in D.items())))

buffers_L, buffers_R = buffers(mz_L, mz_R, max_length=.5)


TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")
_mult = 2
intensity = True
Y_exp = I if intensity else I_d
Y_est = E if intensity else E_d

max_intensity = max(max(Y_est), max(Y_exp))*1.05
plot = figure(plot_width=800*_mult,
              plot_height=400*_mult,
              y_range=(0, max_intensity),
              tools=TOOLS)

plot.xaxis.axis_label = 'mass/charge'
plot.yaxis.axis_label = 'intensity' if intensity else 'intensity / (m over z interval length)'
source = ColumnDataSource(data={'top': Y_exp, 'estimated': Y_est,
                                'left': buffers_L, 'right': buffers_R})

invisible_buffers = plot.quad(left='left',
                              right='right',
                              top='top',
                              bottom=0.0,
                              alpha=0.04,
                              color='black',
                              source=source)
if intensity:
    tolerance = .01
    experimental_bars = plot.vbar(x=masstodon.spectrum.mz,
                                  top=masstodon.spectrum.intensity,
                                  width=tolerance,
                                  color='black',
                                  alpha=.1)

    raw_spectrum = plot.square(x=masstodon.spectrum.mz,
                               y=masstodon.spectrum.intensity,
                               size=10,
                               color='black',
                               alpha=.5)

plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_exp, y1=Y_exp,
             color='black', line_width=3)

plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_est, y1=Y_est,
             color='red', line_width=3)

estimated_bars = plot.quad(left=mz_L,
                           right=mz_R,
                           bottom=0.0,
                           top=Y_est,
                           color='red',
                           fill_alpha=.3)

hover_invisible = HoverTool(renderers=[invisible_buffers],
                            tooltips=[('observed I', "@top{0,0}"),
                                      ('estimated I', "@estimated{0,0}"),
                                      ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
                            mode='vline')

plot.add_tools(hover_invisible)
show(plot)
