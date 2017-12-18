%load_ext autoreload
%autoreload 2

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

sol = masstodon._solutions[0]
sol.node['M0']
sol.node['G0']
sol.node['I0']

sol.node['G0']['intensity']

T_mz = []
T_intensity = []
E_mz_L = []
E_mz_R = []
E_intensity = []
for sol in masstodon._solutions:
    for N in sol:
        if N[0] is 'G':
            E_mz_L.append(sol.node[N]['min_mz'])
            E_mz_R.append(sol.node[N]['max_mz'])
            E_intensity.append(sol.node[N]['intensity'])


# all nodes and edges are tagged with distinct numbers
# masstodon.cz_match.intensities
# masstodon.cz_match.branching_ratio
# masstodon.cz_match.probabilities

G = masstodon._solutions[0]

G.node.data('estimate')



from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool

hover = HoverTool(tooltips=[('intensity', "@top{0,0}"),
                            ('m/z', "[@left{0,0.000}, @right{0,0.000}]")])
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select "
TOOLS += "lasso_select save"
TOOLS = TOOLS.split(' ')
TOOLS.append(hover)

plot = figure(plot_width=300, plot_height=300, tools=TOOLS)
plot.quad(left=np.array([1, 2, 3]),
          right=np.array([1, 2, 3])+.3,
          bottom=0.0,
          top=np.array([1, 2, 3]))
# plot.vbar(left=np.array([1, 2, 3]),
#           right=np.array([1, 2, 3])+.3,
#           top=np.array([1, 2, 3]),
#           color='black')

show(plot)
