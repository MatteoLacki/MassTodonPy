%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from collections import Counter, defaultdict
import numpy as np

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Spectra.Spectrum import Spectrum

substanceP = get_dataset('substanceP')

precursor = {'name': 'substanceP',
             'fasta': substanceP.precursor.fasta,
             'charge': 3}

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
                                               c['estimate'], c['estimate_d'],
                                               c['intensity'], c['intensity_d'])
                                              for (l, r), c in D.items())))

# def get_buffers(mz_L, mz_R, max_length=.45):
#     L = []
#     R = []
#     r_prev = 0  # guards: 0th peak is a phoney
#     for i, (l, r) in enumerate(zip(mz_L, mz_R)):
#         tol = min((l - r_prev)/2, max_length)
#         if i:
#             R.append(r_prev + tol)
#         L.append(l - tol)
#         r_prev = r
#     R.append(r + tol)
#     return L, R
# not symmetric

def get_buffers(mz_L, mz_R, max_length=.45):
    tol = min(max_length, (mz_L[1] - mz_R[0])/2)
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


# source = ColumnDataSource(data=D)
max_intensity = max(max(I_d), max(E_d))
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")
_mult = 2

plot = figure(plot_width=800*_mult,
              plot_height=400*_mult,
              y_range=(0, max_intensity),
              tools=TOOLS)
plot.xaxis.axis_label = 'mass/charge'
plot.yaxis.axis_label = 'intensity'

# tolerance = .03
# experimental_bars = plot.vbar(x=masstodon.spectrum.mz,
#                               top=masstodon.spectrum.intensity,
#                               width=tolerance,
#                               color='black',
#                               alpha=.2)
# experimental_squares = plot.square(x=masstodon.spectrum.mz,
#                                    y=masstodon.spectrum.intensity,
#                                    size=5,
#                                    color='black',
#                                    alpha=.2)

experimental_bars = plot.quad(left=mz_L,
                              right=mz_R,
                              bottom=0.0,
                              top=I_d,
                              line_color='black',
                              fill_alpha=0)

estimated_bars = plot.quad(left=mz_L,
                           right=mz_R,
                           bottom=0.0,
                           top=E_d,
                           color='red',
                           fill_alpha=.3)

invisible_buffers = plot.quad(left=buffers_L,
                              right=buffers_R,
                              top=I_d,
                              bottom=0.0,
                              alpha=0.0)

hover_invisible = HoverTool(renderers=[invisible_buffers],
                            tooltips=[('intensity', "@top{0,0.000}"),
                                      ('m/z', "@left{0,0.000}")],
                            mode='vline')


plot.add_tools(hover_invisible)
show(plot)


#
# plot = figure(plot_width=800*_mult,
#               plot_height=400*_mult,
#               # background_fill_color='black',
#               # border_fill_color='black',
#               y_range=(0, max_intensity),
#               tools=TOOLS)

# hover_quads = HoverTool(renderers = [groups],
#                         tooltips=[('intensity', "@top{0,0}"),
#                                   ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
#                         mode='vline')
# hover_invisible = HoverTool(renderers = [invisible_buffers],
#                             tooltips=[('intensity', "@top{0,0}"),
#                                       ('m/z', "@x{0,0.000}")],
#                             mode='vline')
# plot.add_tools(hover_invisible, hover_bars, hover_squares, hover_quads)

# def get_buffers(mz_L, mz_R, max_buffer_length=.45):
#     r_prev = 0  # guard: 0th peak is a phoney
#     prev_width = 0
#     for i, (l, r) in enumerate(zip(mz_L, mz_R)):
#         buffer_length = min((l- r_prev)/2, max_buffer_length)
#         if i is not 0:
#             yield r_prev - prev_width/2, r_prev + buffer_length # previous right buffer
#         width = r - l
#         yield l - buffer_length, l + width/2  # current left buffer
#         r_prev = r
#         prev_width = width
#     yield r_prev, r_prev + max_buffer_length  # last right buffer
# hover_bars = HoverTool(renderers = [experimental_bars],
#                        tooltips=[('intensity', "@top{0,0}"),
#                                  ('m/z', "[@left{0,0.000}, @right{0,0.000}]")])
#
# hover_bars_e = HoverTool(renderers = [estimated_bars],
#                           tooltips=[('intensity', "@top{0,0}"),
#                                     ('m/z', "[@left{0,0.000},@right{0,0.000}]")])
