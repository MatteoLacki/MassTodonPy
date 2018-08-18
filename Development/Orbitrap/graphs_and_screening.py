%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  defaultdict, namedtuple, Counter
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
from MassTodonPy.Data.Constants      import infinity
from MassTodonPy.models.polynomial   import polynomial
from MassTodonPy.Spectra.lightweight import lightweight_spectrum
from MassTodonPy.Deconvolution.divide_ed_impera import divide_ed_impera

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
min_prob = .8
isotopic_coverage = .99

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
# bc = np.array(list(spec.iter_bc_clusters())) # not needed.
min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads = spec.get_bc_stats()
min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads =\
    [x[good_clusters] for x in (min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads)]

peak_groups = lightweight_spectrum(min_mz, max_mz, total_intensities) # efficient data structure
imperator   = divide_ed_impera(good_mols, peak_groups, min_prob, isotopic_coverage)

imperator.plot()
# add coloring to the edges.



