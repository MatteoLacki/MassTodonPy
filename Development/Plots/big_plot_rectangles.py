%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from bokeh.palettes import viridis
from collections import Counter, defaultdict, namedtuple
from operator import sub
from operator import attrgetter

from Development.Plots.run_masstodon import get_masstodon_results
from MassTodonPy.Plotting.plot_buffers import buffers

masstodon = get_masstodon_results()
sol = masstodon._solutions[0]
mz_digits = masstodon.mz_digits

# TODO:
#   add terminal functions for plotting spectra
#   call ColumnDataSource as many times as there are problems
#   each problem must be represented by a series of basic bricks
#   the data should contain:
#       L, R - constant
#       intensity,                          - constant
#       formula_with_charges, q, g,         - constant
#       parent molecule overall intensity   - constant
#       top, buttom, color                  - adjusted by M.estimate
#
#   substances with larger presence must go at the buttom
#       this should influence the top/buttom values of the bricks
#       I have to calculate first how many intensity there is for each substance
            # that is already stored in sol.node[M]['estimate']
#       or use a priority queue that stores info on the molecules
#
#   TODO:
#       1. write down the methods of the Results class.


# General picture:
# viridis(2)
# ['#440154', '#FDE724']

class ColorGenerator(object):
    def __init__(self):
        self.colors = {}

    def new_color(self):
        """Get the next color in sequence."""
        pass

    def get(self, brick):
        info = (brick.formula, brick.q, brick.g)
        try:
            return self.colors[info]
        except KeyError:
            color = self.new_color()
            self.colors[info] = color
            return color

Cluster = namedtuple('Cluster',
                     'mz_L mz_R mz_left mz_right intensity intensity_d estimate estimate_d')

Brick = namedtuple('Brick', 'cluster top buttom color intensity molecule')

class Results(object):
    """Retrieve results from the MassTodon deconvolution graph."""

    def get_bricks_clusters(self):
        """Generate an unorderd flow of bricks and clusters for the plot."""
        base_width = 10**(-self.masstodon.mz_digits)
        for sol in self.masstodon._solutions:
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
                        cluster = Cluster(mz_L=min_mz,
                                          mz_R=max_mz,
                                          mz_left=min_mz,     # 2 update
                                          mz_right=max_mz,    # 2 update
                                          intensity=intensity,
                                          intensity_d=intensity_d,
                                          estimate=estimate,
                                          estimate_d=estimate_d)
                        yield 'clusters', cluster
                        for I in sol[G]:
                            estimate = sol[G][I]['estimate']
                            for M in sol[I]:
                                if M[0] is 'M':
                                    mol = sol.node[M]['molecule']
                                    yield ('bricks',
                                            Brick(cluster=cluster,
                                                  top=estimate,
                                                  buttom=0.0,
                                                  color='#440154',# 2 update
                                                  intensity=estimate,
                                                  molecule=mol))

    def __init__(self, masstodon):
        self.masstodon = masstodon

        self.bricks = []
        self.clusters = []
        for name, thing in self.get_bricks_clusters():
            self.__dict__[name].append(thing)

        # m/z ratio order
        self.clusters.sort(key=lambda c: c.mz_L)
        # lexicographic order
        self.bricks.sort(key=lambda b: (b.cluster.mz_L, -b.molecule.intensity))

    def adjust_bricks(self, bricks):
        """Adjust the position of the bricks.

        The adjustement is based on their m/z ratio and intensity.
        Also, assign top, buttom, and color to the bricks.
        """
        color = ColorGenerator()
        prev_brick = brick(mz=0.0, intensity=0.0, top=0.0) # guardian
        # while pq not empty:
        #     brick = pq.pop()
        #     brick.color = color.get(brick)
        #     if brick.mz is prev_brick.mz:
        #         brick.bottom = prev_brick.top
        #         brick.top = brick.bottom + brick.intensity
        #     prev_brick = brick
        #     yield brick

    def aggregate(self):
        """Aggregate the results of MassTodon."""
        pass

    def write(self, file_path):
        """Write results to a file.

        Parameters
        ==========
        file_path : str
            Path to the output you want MassTodon to save its results.
            Supported path extensions are .txt and .csv.
        """
        pass

    def plot_in_bokeh(self, file_path=""):
        """Make a bokeh plot.

        Parameters
        ==========
        file_path : str
            Path to the output html file.
            If not provided, MassTodon will use the default temporary location
            used by the Bokeh module.
        """
        pass

results = Results(masstodon)
bricks = results.bricks
clusters = results.clusters

# For making the other plots... probably not so necessary...
mz_L, mz_R = list(zip(*((c.mz_L, c.mz_R) for c in clusters)))
# check out the bokeh flow plotting
mz_lefts, mz_rights = buffers(mz_L, mz_R, max_length=.5)
for c, mz_left, mz_right in zip(clusters, mz_lefts, mz_rights):
    c.mz_left = mz_left
    c.mz_right = mz_right


T = namedtuple('T', 'a b')
t = T(a=[1,2], b='aad')
t.a[0] = 3
t.a
# Find a proper structure for storing mutable thing....
# OK, a list will do. Do it with lists first.
# unless there is a mutable tuple


prev_brick = 0.0
for brick in bricks:
    if brick.mz_L is prev_mz_L:


    prev_mz_L = brick.mz_L



# Old Stuff
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



# class Brick(object):
#     """Container for MassTodon's results."""
#     __slots__ = ('mz_L', 'mz_R', 'mz_left', 'mz_right', 'formula', 'intensity',
#                  'q', 'g', 'mol_intensity', 'mol_name',
#                  'top', 'buttom', 'color')
#
#     def __init__(self, **kwds):
#         """Initialize the brick.
#
#         Pass as many arguments as you know from the '__slots__' tuple,
#         others will be initialized as 'None'.
#         """
#         for k in self.__slots__:
#             setattr(self, k, kwds.get(k, None))
#
#     def __repr__(self):
#         res = "Brick("
#         for k in self.__slots__:
#             res += str(k) + "=" + str(getattr(self, k)) + ", "
#         res = res[0:-2] + ")"
#         return res
#
#     def get_result(self):
#         """Extract the information from the brick."""
#         return self.formula, self.q, self.g, self.intensity
