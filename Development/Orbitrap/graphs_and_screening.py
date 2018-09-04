%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  defaultdict, namedtuple, Counter
import numpy                as      np
import networkx             as      nx
import matplotlib.pyplot    as      plt
from   time                 import  time
import pandas               as      pd
from   math                 import  log10, floor

from MassTodonPy.readers.from_npy               import spectrum_from_npy
from MassTodonPy.Precursor.simple               import precursor
from MassTodonPy.IsotopeCalculator.simple       import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups   import bitonic_clustering
from MassTodonPy.Spectra.simple                 import spectrum
from MassTodonPy.Molecule.simple                import molecule
from MassTodonPy.Formula.Formula                import formula
from MassTodonPy.Molecule.simple                import Molecule
from MassTodonPy.stats.simple_normal_estimators import mean, sd
from MassTodonPy.Data.Constants                 import infinity
from MassTodonPy.models.polynomial              import polynomial
from MassTodonPy.Spectra.lightweight            import lightweight_spectrum
from MassTodonPy.Deconvolution.divide_ed_impera import divide_ed_impera, Imperator, ImperatorMagnus
from MassTodonPy.preprocessing.filters          import filter_subspectra_molecules
from MassTodonPy.Deconvolution.simple           import DeconvolutionProblem
from MassTodonPy.plotters.graphs                import plot_numbered_graph
from MassTodonPy.plotters.spectrum              import plot_spectrum

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.min_mz_diff_clustering()
spec.fit_mz_diff_model()
# spec.plot_mz_diffs()
# spec.plot(clusters='bitonic')
# spec.plot(clusters='min_mz_diff')
spec.bc.stats()



spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()

spec.plot(clusters='bc')
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
isotopic_coverage = .999

good_mols, good_subspectra = filter_subspectra_molecules(subspectra,
                                                         mols,
                                                         std_cnt = 3)
# attention: the sd's will surely change!!! Good! :) Will they? They are not so important.
# The bloody interval widths fully replace this concept.
# to delete
bc = np.array(list(spec.iter_bc_clusters()))
min_mz, max_mz, mean_mz, sds, skewnesses, counts, total_intensities, mz_spreads = spec.get_bc_stats()
ok = min_mz < max_mz
min_mz, max_mz, mean_mz, sds, skewnesses, counts, total_intensities, mz_spreads, bc =\
    [x[ok] for x in (min_mz, max_mz, mean_mz, sds, skewnesses, counts, total_intensities, mz_spreads, bc)]

peak_groups = lightweight_spectrum(min_mz, max_mz, total_intensities) # efficient data structure
t0 = time()
imperator = divide_ed_impera(good_mols, peak_groups, min_prob, isotopic_coverage)
imperator.impera()
fit_time = time() - t0



# imperator.plot()
# imperator.plot_ccs()

ccs    = np.array(imperator.ccs) 
simple = False
cc     = ccs[100] if simple else ccs[np.argmax([len(c) for c in ccs])]







from networkx.linalg.attrmatrix import attr_matrix
# plt_style = 'default'
# plt.style.use(plt_style)

# there is no clear solution to the problem of what should go where.

# this should be given the spectrum, IMHO.
deconvolution_problem

# what if spec contained all these things?
spec[10:400]

dps = []
for cc in ccs:
    dp = DeconvolutionProblem()
    dp.fit(cc, total_intensities, min_mz, max_mz, mean_mz)
    # dp.plot()
    dps.append(dp)

dps = np.array(dps)
dps[10].plot()


fit_to_zeros = False

mz_all   = []
pred_all = []
mz_widths= []
mz_means_all = []
Y_all = []
for cc in ccs:
    mol_columns = np.array([N < 0  for N in cc])
    peak_rows   = np.array([N >= 0 for N in cc])
    X, ordering = attr_matrix(cc, edge_attr='prob')
    X = X[:,mol_columns][peak_rows,:]
    ordering = np.array(ordering)
    peaks    = ordering[ordering >= 0]
    Y = total_intensities[peaks]
    if fit_to_zeros:
        Y = np.concatenate((Y, np.zeros(X.shape[1])))
        x = 1.0 - np.array(X.sum(axis=0)).flatten()
        X = np.concatenate((X,
                            np.diag(x)))
    mz_s = min_mz[peaks]
    mz_e = max_mz[peaks]
    mz_means_all.extend(mean_mz[peaks])
    model = nnls(X, Y)
    beta  = model.coef() 
    pred_all.extend(np.array(np.dot(X[:len(mz_s),], beta)).flatten())
    mz_width = mz_e - mz_s
    mz_widths.extend(mz_e - mz_s)
    mz_all.extend(mz_s)
    Y_all.extend(Y[Y>0])

import matplotlib
# matplotlib.rcParams['figure.figsize'] = 400, 12
plt.bar(mz_all, pred_all, mz_widths,
        align='edge',
        alpha= .5,
        color='grey')
plt.scatter(mz_means_all, Y_all, c= 'red', s=8)
spec.plot(plt_style=plt_style, show=False, peak_color='black')
# plt.savefig('/Users/matteo/Desktop/test.pdf')
plt.show()

# try also not to fit to zero intensities.
# rationale: we alreday filter out some things.
# assume that missigngess is irrelevant now.

betas = []
for cc in ccs:
    Y, X, mz_s, mz_e = get_matrix_representation(cc, total_intensities)
    betas.append(nnls(X, Y))
# this code taks 63.9 ms to solve :D Fuck CVXOPT.





