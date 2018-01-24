%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from bokeh.palettes import viridis
from collections import Counter, defaultdict
from operator import sub

from Development.Plots.run_masstodon import get_masstodon_results
from MassTodonPy.Plotting.plot_buffers import buffers

# TODO:
#   parallelize the making of the plot if it takes a lot of time
#       check ubiquitin example
#   [safari does not read the plots.] SOLVED
#   add terminal functions for plotting spectra

masstodon = get_masstodon_results()
sol = masstodon._solutions[0]
mz_digits = masstodon.mz_digits

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
            if intensity > 0:
                intensity_d = intensity/(max_mz - min_mz)
                estimate = sol.node[G]['estimate']
                estimate_d = estimate/(max_mz - min_mz)
                out_counter = Counter({'intensity': intensity,
                                       'intensity_d': intensity_d,
                                       'estimate': estimate,
                                       'estimate_d': estimate_d})
                out_formulas = {}
                for I in sol[G]:
                    GM_estimate = sol[G][I]['estimate']
                    for M in sol[I]:
                        if M[0] is 'M':
                            mol = sol.node[M]['molecule']
                            mol = (str(mol.formula), mol.q, mol.g)
                            out_formulas[mol] = GM_estimate
                yield (min_mz, max_mz), (out_counter, out_formulas)


def get_data_for_one_solution(sol, mz_digits):
    precise = defaultdict(dict)
    aggregated = defaultdict(Counter)
    for interval, (data, formulas) in get_info_on_solution(sol, mz_digits):
        aggregated[interval] += data
        precise[interval].update(formulas)
    molecule_tags = set([])
    for interval in precise:
        for mol, intensity in precise[interval].items():
            molecule_tags.add(mol)
    molecule_tags = list(molecule_tags)
    data = defaultdict(list)
    for L, R in precise:
        for mol in molecule_tags:
            data[mol].append(precise[(L, R)].get(mol, 0.0))
        data['mz'].append((L + R)/2.0)
        data['width'].append(R - L)
    data_str = {}
    for k, v in data.items():
        data_str[str(k)] = v
    return aggregated, molecule_tags, data_str

output_file("big_plot.html", mode='inline')
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")
_mult = 1
intensity = True
plot = figure(plot_width=800*_mult,
              plot_height=400*_mult,
              tools=TOOLS)
plot.y_range.start = 0
plot.xaxis.axis_label = 'mass/charge'
plot.yaxis.axis_label = 'intensity'
tolerance = 0.01

experimental_bars = plot.vbar(x=masstodon.spectrum.mz,
                              top=masstodon.spectrum.intensity,
                              width=tolerance,
                              color='black',
                              alpha=.1)

raw_spectrum = plot.square(x=masstodon.spectrum.mz,
                           y=masstodon.spectrum.intensity,
                           size=5,
                           color='black',
                           alpha=.5)


sol = masstodon._solutions[0]
aggregated = defaultdict(Counter)
data = []
mol_tags = []
for sol in masstodon._solutions:
    agg, _mol_tags, _data = get_data_for_one_solution(sol, mz_digits)
    for interval, counter in agg.items():
        aggregated[interval] += counter
    data.append(_data)
    mol_tags.append(_mol_tags)

mz_L, mz_R, I, I_d, E, E_d = list(zip(*sorted((
    l, r,
    c['intensity'],
    c['intensity_d'],
    c['estimate'],
    c['estimate_d']) for (l, r), c in aggregated.items())))

buffers_L, buffers_R = buffers(mz_L, mz_R, max_length=.5)
buffer_map = {}
for L, R, bL, bR in zip(mz_L, mz_R, buffers_L, buffers_R):
    buffer_map[(L + R)/2.0] = (bL, bR)

def get_colors(molecule_tags):
    return viridis(len(molecule_tags))

_mol_tags, _data = mol_tags[2], data[2]
for i, (_mol_tags, _data) in enumerate(zip(mol_tags, data)):
    try:
        colors = get_colors(_mol_tags)
        _data['buffers'] = [ -sub(*buffer_map[mz]) for mz in _data['mz'] ]
        source = ColumnDataSource(_data)
        _tags = [str(m) for m in _mol_tags]
        plot.vbar_stack(stackers=_tags, x='mz', width='buffers',
                        color=colors, source=source, alpha=0.04)
        plot.vbar_stack(stackers=_tags, x='mz', width='width',
                        color=colors, source=source)

        # hover_invisible = HoverTool(renderers=[invisible_buffers],
        #                             tooltips=[('observed I', "@top{0,0}"),
        #                                       ('estimated I', "@estimated{0,0}"),
        #                                       ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
        #                             mode='vline')
        #
        # plot.add_tools(hover_invisible)

    except:
        print(i)

plot.segment(x0=buffers_L, x1=buffers_R, y0=I, y1=I,
             color='black', line_width=3)

plot.segment(x0=buffers_L, x1=buffers_R, y0=E, y1=E,
             color='red', line_width=3)



show(plot)






# hover_invisible = HoverTool(renderers=[invisible_buffers],
#                             tooltips=[('observed I', "@top{0,0}"),
#                                       ('estimated I', "@estimated{0,0}"),
#                                       ('m/z', "[@left{0,0.000}, @right{0,0.000}]")],
#                             mode='vline')
#
# plot.add_tools(hover_invisible)




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
