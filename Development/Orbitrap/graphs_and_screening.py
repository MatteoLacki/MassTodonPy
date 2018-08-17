%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  Counter
import numpy                as      np
import matplotlib.pyplot    as      plt
from   time                 import  time
import intervaltree         as      iTree
import networkx             as      nx
from   networkx             import  connected_component_subgraphs
import pandas               as      pd
from   math                 import  log10, floor

from MassTodonPy.readers.from_npy             import spectrum_from_npy
from MassTodonPy.Precursor.simple             import precursor
from MassTodonPy.IsotopeCalculator.simple     import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups import bitonic_clustering
from MassTodonPy.Spectra.simple               import spectrum
from MassTodonPy.Molecule.simple              import molecule
from MassTodonPy.Formula.Formula              import formula
from MassTodonPy.Molecule.simple              import Molecule
from MassTodonPy.stats.simple_normal_estimators   import mean,\
                                                         sd
from MassTodonPy.Data.Constants     import infinity
from MassTodonPy.models.polynomial  import polynomial



# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()
# TODO: replace this later on: a more complex model has to be fitted.
# and the estimation of the Standard Deviations of groups has to be 
# dependent upon more data points all across the spectrum.
spec.fit_sd_mz_model()

subspectra = list(spec.iter_mdc_subspectra())
mz_digits  = spec.get_mz_digits()
iso_calc   = isotope_calculator(digits = mz_digits)

# generating formulas
fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name="", iso_calc=iso_calc)
mols   = np.array(list(prec.molecules()))

# Screen stupid formulas, and then proceed with standard algorithm for 
# finding connected components of the neighbourhood graph.
emp_tree = iTree.IntervalTree()
for subspec in subspectra:
    s, e = subspec.interval
    emp_tree[s:e] = subspec

good_mols = []
good_subspectra = set([])
good_subspectra_cnt = Counter()
for mol in mols:
    s, e = mol.interval(std_cnt = 3)
    touched_spectra = emp_tree[s:e]
    if touched_spectra:
        good_mols.append(mol)
        good_subspectra |= touched_spectra



# attention: the sd's will surely change!!! Good! :)
#to delete
bc = np.array(list(spec.iter_bc_clusters()))
# when to eliminate them sinlgy peaked clusters???

min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads = spec.get_bc_stats()
good_clusters = min_mz < max_mz

# instead of filtering out these 3 peaks, we could modify them:
#   add +- estimated delta m/z to left and right ends.
# hopefully, this will not upset the differences, but who knows.
min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads =\
    [x[good_clusters] for x in (min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads)]



from bisect import bisect_left, bisect_right, §


class lightweight_spectrum(object):
    def __init__(self, min_mz, max_mz, total_intensities):
        self.clust_no = len(min_mz)
        self._spec = np.zeros(dtype = float,
                              shape = (self.clust_no*2,))
        # real spectrum
        self._spec[0::2] = min_mz
        self._spec[1::2] = max_mz
        self.group_intensity = total_intensities

    def __getitem__(self, key):
        i = bisect(self._spec, key)
        i_div, i_mod = divmod(i, 2)
        print(i_mod)
        return i_div if i_mod else -1

    def __repr__(self):
        return self._spec.__repr__()


L = lightweight_spectrum(min_mz, max_mz, total_intensities)
L = lightweight_spectrum([0,5,8], 
                         [1,6,9],
                         total_intensities)
L[10]
L[-2]
L[0]
L[.5]
L[1]
L[1.5]
L[4]
L[5]
L[5.5]
L[6]

lspec._spec.shape
lspec._spec

last_e = lspec._spec[-1,0]

value = 3000
i = bisect(lspec._spec[:,0], value)
i_div_2, i_mod_2 = divmod(i, 2)



i -= 1
lspec._spec[i]
if i_mod_2 == 0:
    # we are in [E_{i_div_2},S_{i_div_2 + 1})
    # we don't need to do anything here.
    return None
else:
    return i-1
lspec._spec[(i-1):(i+1),:]


i = bisect(lspec._spec[:,0], 3000)
outer_e = lspec._spec[i-1,0]
outer_e
bisect(lspec._spec[:,0], outer_e)
i

bisect(lspec._spec[:,0], (outer_e + lspec._spec[i-2,0])/2.0)

if i == len(lspec._spec):
    # we are beyond the 

lspec._spec[i,0]



key = 200
i = bisect(lspec._spec[:,0], key)

lspec._spec[i]





bisect_left(clusters_info[:,0], 10)
bisect_right(clusters_info[:,0], 10)

clusters_info[0:4,0]
bisect_left(clusters_info[0:4,0], 150.69678862)
bisect_right(clusters_info[0:4,0], 150.69678862)

bisect_left(clusters_info[0:4,0], 150.69678863)
bisect_right(clusters_info[0:4,0], 150.69678863)

# understanding bisect_left and bisect right
x = np.array([0,1, 4,5, 8,9])

bisect_right(x, -1)
bisect_right(x, 0)
bisect_right(x, 0.5)
bisect_right(x, 1)
bisect_right(x, 8.5)
bisect_right(x, 9)
# nieparzyste: wpadamy w przedział.
# parzyste: nie wpadamy w przedział, chyba że trafiliśmy w lewy kraniec.




bisect_right(x, 0)
bisect_right(x, 1.5)
x
x[bisect_right(x, 4)]

bisect(x, .5)
bisect_right(x, -1)





%%timeit
mol_G = nx.Graph()
# constriuction takes 4.5 micro-sec: not bad.


bc_clusters = np.array(list(spec.iter_bc_clusters()))


for mol in good_mols:
    # do I actually need this graph? 
    mol_G = nx.Graph()













