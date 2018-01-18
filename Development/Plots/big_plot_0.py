%load_ext autoreload
%autoreload 2

from collections import Counter, defaultdict

from Development.Plots.run_masstodon import get_masstodon_results
from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool

from MassTodonPy.Plotting.plot_buffers import buffers

masstodon = get_masstodon_results()
sol = masstodon._solutions[0]

def get_info_on_solution(sol, mz_digits):
    base_width = 10**(-mz_digits)
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


E = defaultdict(list)
D = defaultdict(Counter)
for interval, (data, formulas) in get_info_on_solution(sol, mz_digits):
    D[interval] += data
    E[interval].extend(formulas)

mz_L, mz_R, I, I_d, E, E_d = list(zip(*sorted((l, r,
                                               c['intensity'],
                                               c['intensity_d'],
                                               c['estimate'],
                                               c['estimate_d'])
                                              for (l, r), c in D.items())))
buffers_L, buffers_R = buffers(mz_L, mz_R, max_length=.5)




# Big Idea:
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")
_mult = 2
intensity = True
plot = figure(plot_width=800*_mult,
              plot_height=400*_mult,
              tools=TOOLS)
plot.y_range.start = 0
plot.xaxis.axis_label = 'mass/charge'
plot.yaxis.axis_label = 'intensity'

for sol in masstodon._solutions:
    molecule_tags, source = get_data_for_one_solution(sol, mz_digits)
    # source must contain all fields from molecule_tags !!!
    colors = get_colors(molecule_tags)
    plot.vbar_stack(molecule_tags,
                    x='mz',
                    width='widths',
                    color=colors,
                    source=source)

plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_exp, y1=Y_exp,
             color='black', line_width=3)

plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_est, y1=Y_est,
             color='red', line_width=3)

hover_invisible = HoverTool(renderers=[invisible_buffers],
                            tooltips=[('observed I', "@top{0,0}"),
                                      ('estimated I', "@estimated{0,0}"),
                                      ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
                            mode='vline')

plot.add_tools(hover_invisible)
show(plot)











# source = ColumnDataSource(data={'top': Y_exp, 'estimated': Y_est,
#                                 'left': buffers_L, 'right': buffers_R})
#
# invisible_buffers = plot.quad(left='left',
#                               right='right',
#                               top='top',
#                               bottom=0.0,
#                               alpha=0.04,
#                               color='black',
#                               source=source)
# if intensity:
#     tolerance = .01
#     experimental_bars = plot.vbar(x=masstodon.spectrum.mz,
#                                   top=masstodon.spectrum.intensity,
#                                   width=tolerance,
#                                   color='black',
#                                   alpha=.1)
#
#     raw_spectrum = plot.square(x=masstodon.spectrum.mz,
#                                y=masstodon.spectrum.intensity,
#                                size=10,
#                                color='black',
#                                alpha=.5)
#
# plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_exp, y1=Y_exp,
#              color='black', line_width=3)
#
# plot.segment(x0=buffers_L, x1=buffers_R, y0=Y_est, y1=Y_est,
#              color='red', line_width=3)
#
# estimated_bars = plot.quad(left=mz_L,
#                            right=mz_R,
#                            bottom=0.0,
#                            top=Y_est,
#                            color='red',
#                            fill_alpha=.3)
#
# hover_invisible = HoverTool(renderers=[invisible_buffers],
#                             tooltips=[('observed I', "@top{0,0}"),
#                                       ('estimated I', "@estimated{0,0}"),
#                                       ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
#                             mode='vline')
#
# plot.add_tools(hover_invisible)
# show(plot)
