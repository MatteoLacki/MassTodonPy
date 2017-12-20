%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from collections import Counter, defaultdict
import numpy as np

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
# from MassTodonPy.Misc.plot_buffers import get_buffers
from MassTodonPy.Misc.sorting import sort_by_first
from MassTodonPy.Spectra.Spectrum import Spectrum

substanceP = get_dataset('substanceP')
spectrum = substanceP.spectrum

spectrum.plot()

atoms = spectrum.atoms
masses = spectrum.masses
buffers_L, buffers_R = get_buffers(atoms, atoms, max_length=2)
plot_width = 1000
plot_height = 600
bar_width = .01

TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(' ')
plot = figure(tools=TOOLS,
              width=plot_width,
              height=plot_height)

_store_names = spectrum._store_names

source = ColumnDataSource(data={'top': masses,
                                'left': buffers_L,
                                'right': buffers_R,
                                'atoms': atoms})

invisible_buffers = plot.quad(left='left',
                              right='right',
                              top='top',
                              bottom=0.0,
                              alpha=0.0,
                              color='black',
                              source=source)
raw_spectrum = plot.vbar(x=atoms,
                         top=masses,
                         width=bar_width,
                         color='black')
hover = HoverTool(renderers=[invisible_buffers],
                  tooltips=[(_store_names[1], "@top{0,0}"),
                            (_store_names[0], "@atoms{0,0.000}")],
                  mode='vline')
plot.add_tools(hover)
show(plot)





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

def base_data(solutions, mz_digits):
    base_width = 10**(-mz_digits)
    for sol in solutions:
        for N in sol:
            if N[0] is 'G':
                try:
                    min_mz = sol.node[N]['min_mz'] - base_width/2
                    max_mz = sol.node[N]['max_mz'] + base_width/2
                except KeyError:
                    mz = sol.node[N]['mz']
                    min_mz = mz - base_width/2
                    max_mz = mz + base_width/2
                    # pass
                intensity = sol.node[N]['intensity']
                if intensity > 0 or estimate > 0:
                    intensity_d = intensity/(max_mz - min_mz)
                    estimate = sol.node[N]['estimate']
                    estimate_d = estimate/(max_mz - min_mz)
                    yield (min_mz, max_mz), Counter({'intensity': intensity,
                                                     'intensity_d': intensity_d,
                                                     'estimate': estimate,
                                                     'estimate_d': estimate_d})

D = defaultdict(Counter)
for interval, data in base_data(masstodon._solutions, masstodon.mz_digits):
    D[interval] += data

mz_L, mz_R, I, I_d, E, E_d = list(zip(*sorted((l, r,
                                               c['intensity'], c['intensity_d'],
                                               c['estimate'], c['estimate_d'])
                                              for (l, r), c in D.items())))

buffers_L, buffers_R = get_buffers(mz_L, mz_R, max_length=.5)


buffers_L, buffers_R = list(zip(*get_buffers(mz_L, mz_R)))


def get_buffers(mz_L, mz_R, max_length=.5):
    r_prev = -infinity
    l = mz_R[0]
    r = mz_R[0]
    l_next = mz_L[1]
    tol = min(max_length, (l-r_prev)/2, (l_next - r)/2)
    L = [mz_L[0]-tol]
    R = [mz_R[0]+tol]
    for i in range(1,len(mz_L)-1):
        r_prev = mz_R[i-1]
        l = mz_L[i]
        r = mz_R[i]
        l_next = mz_L[i+1]
        tol = min(max_length, (l-r_prev)/2, (l_next-r)/2)
        L.append(l-tol)
        R.append(r+tol)
    tol = min(max_length, (mz_L[-1] - mz_R[-2])/2)
    L.append(mz_L[-1] - tol)
    R.append(mz_R[-1] + tol)
    return L, R

buffers_L, buffers_R = get_buffers(mz_L, mz_R)
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
