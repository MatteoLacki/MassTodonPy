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
from MassTodonPy.Molecule.simple              import Molecule
from MassTodonPy.stats.simple_normal_estimators   import mean,\
                                                         sd

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()
spec.fit_sd_mz_model()
spec.sd_mz_model.plot()

subspectra = list(spec.iter_mdc_subspectra())

spec[1126.5:1126.8].plot()
spec[1756:1766].plot()

means, sds, counts, total_intensities, mz_spreads = spec.get_bc_stats()

DD = pd.DataFrame(dict(mean=means, sd=sds, intensity=total_intensities, cnt=counts, spread=mz_spreads))
DD.to_csv(path_or_buf = '/Users/matteo/Projects/MassTodonPy/data/data_4_sd.csv', index=False)


# generating formulas
fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name = "")
mols   = list(prec.molecules())

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
    s, e    = mol.interval(std_cnt = 3)
    touched_spectra = emp_tree[s:e]
    if touched_spectra:
        good_mols.append(mol)
        good_subspectra |= touched_spectra
# the screening lasts: 4.12 sec
# this naive filtering managed to cut the problem by more than half.
len(good_mols)
len(good_subspectra)

# test if its worth it to use the prerviously established clusters.

from MassTodonPy.Deconvolution.Deconvolve import build_deconvolution_graph


# simplify the isotopic distributions: go back to merging isotopologues.

# is it worth it?
DG = build_deconvolution_graph(good_mols, spec)



