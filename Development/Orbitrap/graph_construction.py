%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  Counter
import numpy                as      np
import matplotlib.pyplot    as      plt
from   time                 import  time
import intervaltree         as      iTree
import networkx             as      nx
from networkx               import  connected_component_subgraphs

from MassTodonPy.readers.from_npy             import spectrum_from_npy
from MassTodonPy.Precursor.simple             import precursor
from MassTodonPy.IsotopeCalculator.simple     import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups import bitonic_clustering
from MassTodonPy.Spectra.simple               import spectrum
from MassTodonPy.Molecule.simple              import molecule
from MassTodonPy.Formula.Formula              import formula


# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()
subspectra = list(spec.iter_mdc_subspectra())

# generating formulas
fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name = "")
mols   = list(prec.molecules())

# mol_tree = iTree.IntervalTree()
# for mol in mols:
#     s, e          = mol.interval 
#     mol_tree[s:e] = mol
# 10 sec: long long

emp_tree = iTree.IntervalTree()
for subspec in subspectra:
    s, e = subspec.interval
    emp_tree[s:e] = subspec
# 25.3 ms ± 1.25 ms per loop (mean ± std. dev. of 7 runs, 100 loops each)

# this seems a lot quicker than doing the ccs
dd = {}
for mol in mols:
    s, e    = mol.interval
    dd[mol] = emp_tree[s:e] 
# but this can be done much faster by direct iteration from left to right!

# alternative to the above: construct a neighbourhood graph.
def build_graph_with_interval_tree(iTree, mols):
    G = nx.Graph()
    for M in mols:
        s, e = M.interval
        G_M = iTree[s:e]
        if G_M:
            G.add_node(M, type = 'M')
        for S in G_M:
            G.add_node(S, type = 'S')
            G.add_edge(S, M)
    return G

G = build_graph_with_interval_tree(emp_tree, mols)
# 8.2 secs

cs = list(connected_component_subgraphs(G))
Counter(len(c) for c in cs)

spec[759.3382795685296:760.5482233729739].plot()
Counter(len(v) for k, v in dd.items())
# Counter({0: 33074, 1: 29335, 2: 2077, 3: 95, 4: 5}): hence, it seems that mostly things do not compete.
# have to check, if the bloody intervals are calculated properly.

# another approach: build a graph based on the established bc groups.
bc_tree = iTree.IntervalTree()
for subspec in subspectra:
    for bc_group in subspec.iter_bc_subspectra():
        s, e = bc_group.interval
        if s < e: # what to do otherwise: most likely there is an error.
            bc_tree[s:e] = bc_group
# 13.8 sec
G = build_graph_with_interval_tree(bc_tree, mols)
# 64 secs
cs = list(connected_component_subgraphs(G))
#  96 secs
Counter(len(c) for c in cs)
# the analysis looks very similar, like the previous one.



