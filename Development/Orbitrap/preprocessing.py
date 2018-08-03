%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import numpy as np
import matplotlib.pyplot as plt

from MassTodonPy.plotters.spectrum              import plot_spectrum
from MassTodonPy.readers.from_npy               import spectrum_from_npy
from MassTodonPy.Spectra.peak_clustering        import mz_bitonic
from MassTodonPy.Spectra.peak_clustering        import iter_cluster_ends
from MassTodonPy.Spectra.orbitrap_peak_groups   import bitonic_clustering
from MassTodonPy.models.polynomial              import polynomial
from MassTodonPy.models.spline                  import spline
from MassTodonPy.stats.simple_normal_estimators import mean, sd
from MassTodonPy.plotters.spectrum              import plot_peak_group


data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)
bc            = bitonic_clustering(mz,
                                   intensity, 
                                   min_mz_diff  = .15,
                                   abs_perc_dev = .2)
bc.plot_mz_diffs()
# bc.fit_mz_diffs(model = polynomial)
# bc.plot_mz_diffs()

# # the fishy m/z values are here.
# fishy_mz = bc.mz_diff_model.x[ ~bc.mz_diff_model.is_signal ]
# fishy_clusters = np.unique(bc.clusters[np.isin(bc.mz, fishy_mz)])
# fc = fishy_clusters[3]

# def plot_fishy(fc):
#     mz_fishy        = mz[bc.clusters == fc]
#     intensity_fishy = intensity[bc.clusters == fc]

#     plt.vlines(mz_fishy, [0], intensity_fishy, colors='white')
#     plt.show()

# plot_fishy(fishy_clusters[7])
# # ergo: not all fishy clusters are so fishy. Leave it for now.


### The estimation of thłe Gaussian parameters.
# INPUT:    binned histogram data
# OUTPUT:   estimated mean and standard deviation of the normal distribution


# this does not work at all: all the intensities are almost the same...

# Code below look as some methods for an OrbitrapSpectrum class.
groups = list(bc.iter_clusters())
mz_g, intensity_g = groups[1500]

# the overall plot:
means = np.array([mean(mz, i) for mz, i in bc.iter_clusters()])
sds   = np.array([sd(mz, i) for mz, i in bc.iter_clusters()])

len(mz) / len(means)


OK = np.ones(shape = means.shape) * .2
OK[[67, 137, 138, 207, 238, 290, 314, 357]] = 5.0
plt.style.use('dark_background')
plt.scatter(means, sds, s = OK, c = OK==.2)
plt.show()

s = spline(means, sds)
s.plot()
p = polynomial(means, sds)


# cluster_no = 0
# groups_iterator = bc.iter_clusters()
# nrows = 3
# ncols = 6
# for i in range(nrows * ncols):
#     plt.subplot(nrows, ncols, i+1)
#     try:
#         plot_peak_group(*next(groups_iterator), show=False)
#     except Exception:
#         pass
#     plt.title("Cluster No {}".format(cluster_no))
#     cluster_no += 1
# plt.show()
from MassTodonPy.Measure.Measure import Measure

m = Measure(atoms=mz, masses=intensity)


