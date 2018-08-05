%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from time import time


from MassTodonPy.readers.from_npy    import spectrum_from_npy
from MassTodonPy.Precursor.Precursor import precursor
from MassTodonPy.IsotopeCalculator.simple           import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups       import bitonic_clustering
from MassTodonPy.Spectra.simple                     import spectrum

data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)
spec = spectrum(mz, intensity)
spec.plot()
spec.bitonic_clustering()
spec.plot()
spec.min_mz_diff_clustering()
spec.plot()

spec.mz_lefts_mz_diffs_in_clusters()
spec.fit_mz_diff_model()
spec.mz_diff_model(10)
spec.plot_mz_diffs()

# deduplication should be in the spectrum
from MassTodonPy.Spectra.orbitrap.peak_clustering import ClusteringAlgorithm


fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name = "shit")
mols   = list(prec.molecules())


spectrum(mz, intensity)

len(mols)
# exchange the molecule's isotopic generator for simpler.
mol = mols[0]
mol.mean_mz
mol.monoisotopic_mz
f = mol.formula
# Task: enhance the calculation
# t0 = time()
# isotopologues = [m.isotopologues() for m in mols]
# t1 = time()
# 89 secs.
# this is actually not so bad,
# the question is, do we need to do it?
# but how much can we save on using mean and sd?
# OK, have to have it done in Numpy to have a speed-up
iso_calc = isotope_calculator()
iso_calc.monoisotopic_mz(mol.formula, q=charge)
iso_calc.heaviest_mz(mol.formula, q=charge) - iso_calc.lightiest_mz(mol.formula, q=charge)

%%timeit
iso_calc.most_probable_mz(mol.formula, q=charge)


%%timeit
mean = iso_calc.mean_mz(mol.formula, q=charge)
sd   = iso_calc.sd(mol.formula,   q=charge)


mean - 3*sd < env.head_mz() / charge
mean + 3*sd > env.tail_mz() / charge

env = iso_calc(mol.formula)
env.plot()
env.tail_mz()
env.head_mz()
(env.tail_mz() - env.head_mz())/charge
env.max_peak()

# get a method that defines the maximal space for the given
# set of atoms in the precursor

max(max(np.diff(iso_calc._masses[e])) for e in prec.formula)

# how much can we gain by subdivision of spectrum based on m/z ?
# we can actually try to do bitonic_clustering after this obvious preprocessing.


np.diff(iso_calc._masses['C'])
np.diff(iso_calc._masses['H'])
np.diff(iso_calc._masses['N'])
np.diff(iso_calc._masses['S'])
np.diff(iso_calc._masses['O'])
# so, a good value is 1.1 Thomson. I wouldn't believe that it can hurt.
max_diff = 1.1


list(max_mz_diff_iterator(mz, max_diff))[-1]
len(list(iter_cluster_ends(max_mz_diff_iterator(mz, max_diff)))

for mz_l, intensity_l in iter_clusters(mz,
                                       intensity,
                                       max_mz_diff_iterator(mz, 1.5)):

subspectra = [ 
    for spec in iter_clusters(mz, intensity, max_mz_diff_iterator(mz, max_diff))]


%%timeit
env = iso_calc(mol.formula, _memoize=False)




%%timeit
means = [iso_calc.mean(m.formula) for m in mols]




iso_calc.get_envelope('C200H400')
iso_calc('C200H400')



atoms_in_all_formulas = set(prec.formula._storage.keys())



mol.formula
A = iso_calc(mol.formula)
masses, probs = A.atoms, A.masses


plot_spectrum(masses, probs, show=False)
plt.scatter(mol.mean_mz, 0, c= 'red')
plt.show()

