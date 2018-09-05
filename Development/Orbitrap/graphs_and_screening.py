%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from   collections          import  defaultdict, namedtuple, Counter
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
from MassTodonPy.Deconvolution.divide_ed_impera import divide_ed_impera
from MassTodonPy.preprocessing.filters          import filter_subspectra_molecules
from MassTodonPy.Deconvolution.simple           import DeconvolutionProblem
from MassTodonPy.plotters.graphs                import plot_numbered_graph
from MassTodonPy.plotters.spectrum              import plot_spectrum

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering(model_diff=None, model_sd=None)
spec.min_mz_diff_clustering()

# spec.plot_mz_diffs()
# spec.plot(clusters='bitonic')
# spec.plot(clusters='min_mz_diff')
# spec.bc.plot_sd()
# spec.plot(clusters='bitonic')
# spec.plot(clusters='min_mz_diff')

subspectra = list(spec.iter_min_mz_diff_subspectra())
mz_digits  = spec.bc.get_smallest_diff_digits()
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
t0 = time()
imperator = divide_ed_impera(good_mols, spec.bc, min_prob, isotopic_coverage)
fit_time = time() - t0
# imperator.plot()
# imperator.plot_ccs()

imperator.set_estimated_intensities()




# what is this code below????
    # a global view on spectrum after fitting
        # where should it be? which class it belongs to?
            # spectrum?


# it seems we have to now divide the fittings according to some criterion 
# into shit, soso, hmmm, ok, good, great.
from MassTodonPy.models.nnls import nnls
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
plt_style = 'fast'
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





