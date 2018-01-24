%load_ext autoreload
%autoreload 2

from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool
from bokeh.palettes import viridis, Colorblind
from collections import Counter, defaultdict, namedtuple
from itertools import cycle
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

class Cluster(object):
    __slots__ = ('mz_L', 'mz_R', 'mz_left', 'mz_right', 'intensity',
                 'intensity_d', 'estimate', 'estimate_d')

    def __init__(self, **kwds):
        for s in self.__slots__:
            setattr(self, s, kwds.get(s, None))

    def __repr__(self):
        res = []
        for s in self.__slots__:
            v = getattr(self, s)
            res.append( "{0}={1}".format(s,v))
        return "Cluster(" + ", ".join(res) + ")"


class Brick(object):
    __slots__ = ('cluster', 'top', 'buttom', 'color', 'intensity', 'molecule')

    def __init__(self, **kwds):
        for s in self.__slots__:
            setattr(self, s, kwds.get(s, None))

    def __repr__(self):
        res = []
        for s in self.__slots__:
            v = getattr(self, s)
            res.append( "{0}={1}".format(s,v))
        return "Brick(" + ", ".join(res) + ")"


class ColorGenerator(object):
    """Generator of colors with memoization."""

    def __init__(self, max_colors=8):
        self.color = cycle(Colorblind[max_colors])
        self.colors = {}

    def __call__(self, molecule):
        """Get the next color in sequence."""
        try:
            return self.colors[molecule]
        except KeyError:
            self.colors[molecule] = next(self.color)
            return self.colors[molecule]


class Results(object):
    """Retrieve results from the MassTodon deconvolution graph."""

    def __init__(self, masstodon, max_buffer_len=0.5, max_colors=8):
        self.masstodon = masstodon
        self.bricks = []
        self.clusters = []
        for name, thing in self.get_bricks_clusters():
            self.__dict__[name].append(thing)
        self.clusters.sort(key=lambda c: c.mz_L) # order by m/z ratios
        # adding left and right buffers needed for plot
        mz_L, mz_R = list(zip(*((c.mz_L, c.mz_R) for c in self.clusters))) #TODO this is an overkill: smooth out the defintion of the buffers function.
        mz_lefts, mz_rights = buffers(mz_L, mz_R, max_length=max_buffer_len)
        for c, mz_left, mz_right in zip(self.clusters, mz_lefts, mz_rights):
            c.mz_left = mz_left
            c.mz_right = mz_right

        self.color = ColorGenerator(max_colors)
        # lexicographically order the bricks
        self.bricks.sort(key=lambda b: (b.cluster.mz_L, -b.molecule.intensity))
        # adjust brick heights
        prev_b = Brick(top=0.0, buttom=0.0, cluster=Cluster(mz_L=0.0)) # guardian
        for b in self.bricks:
            if b.cluster.mz_L is prev_b.cluster.mz_L:
                b.buttom = prev_b.top
                b.top = prev_b.top + b.intensity
            b.color = self.color(b.molecule)
            prev_b = b


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

    def adjust_bricks(self, bricks):
        """Adjust the position of the bricks.

        The adjustement is based on their m/z ratio and intensity.
        Also, assign top, buttom, and color to the bricks.
        """
        pass

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

    def plot_in_bokeh(self, file_path="", mode="inline"):
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


# Making the Plot

file_path = "big_plot.html"
mode = 'inline'
output_file(file_path, mode=mode)
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select save".split(" ")
_mult = 1
plot = figure(plot_width=800*_mult,
              plot_height=400*_mult,
              tools=TOOLS)
plot.y_range.start = 0
plot.xaxis.axis_label = 'mass/charge'
plot.yaxis.axis_label = 'intensity'
tolerance = 0.01
mz_representation = "@mz{0." + "".join(["0"] * mz_digits) + "}"

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

# Plotting the bricks
lists = zip(*((b.cluster.mz_L,
               b.cluster.mz_R,
               b.cluster.mz_left,
               b.cluster.mz_right,
               (b.cluster.mz_L + b.cluster.mz_R)/2.0,
               b.buttom,
               b.top,
               b.color,
               b.intensity,
               b.molecule.name,
               b.molecule.source,
               str(b.molecule.formula),
               b.molecule.q,
               b.molecule.g,
               b.molecule.intensity)
              for b in bricks))

data = dict(zip(('mz_L',
                 'mz_R',
                 'mz_left',
                 'mz_right',
                 'mz',
                 'buttom',
                 'top',
                 'color',
                 'intensity',
                 'molname',
                 'molsource',
                 'molformula',
                 'molq',
                 'molg',
                 'molintensity'),
                lists))

source = ColumnDataSource(data)

fat_rectangles = plot.quad(top='top', bottom='buttom',
                           left='mz_left', right='mz_right',
                           color='color', source=source,
                           alpha=.5)

slim_rectangles = plot.quad(top='top', bottom='buttom',
                           left='mz_L', right='mz_R',
                           color='color', source=source)

hover_fat = HoverTool(renderers=[fat_rectangles],
                      tooltips=[('peak intensity', "@intensity{0,0}"),
                                ('m/z interval', mz_representation),
                                ('name', '@molname'),
                                ('formula','@molformula'),
                                ('q', '@molq'),
                                ('g', '@molg'),
                                ('total intensity', "@molintensity{0,0}")])

plot.add_tools(hover_fat)

# Adding aggregate information

cluster_lists = zip(*((c.mz_left,
                       c.mz_right,
                       (c.mz_left + c.mz_right)/2.0,
                       c.intensity,
                       c.estimate,
                       abs(c.intensity-c.estimate))
                      for c in clusters))

data_clusters = dict(zip(('mz_left',
                          'mz_right',
                          'mz',
                          'intensity',
                          'estimate',
                          'abserror'),
                         cluster_lists))

source_clusters = ColumnDataSource(data_clusters)

cluster_intensities = plot.segment(x0='mz_left',
                                   x1='mz_right',
                                   y0='intensity',
                                   y1='intensity',
                                   color='red',
                                   line_width=3, source=source_clusters)

hover_clusters = HoverTool(renderers=[cluster_intensities],
                           tooltips=[('total intensity', "@intensity{0,0}"),
                                     ('total estimate', "@estimate{0,0}"),
                                     ('absolute error', "@abserror{0,0}"),
                                     ('m/z', mz_representation)])

plot.add_tools(hover_clusters)

show(plot)
