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


data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)
bc            = bitonic_clustering(mz,
                                   intensity, 
                                   in_mz_diff   = .15,
                                   abs_perc_dev = .2)

bc.fit_mz_diffs()
bc.plot_mz_diffs()
bc.fit_mz_diffs(model = polynomial)
bc.plot_mz_diffs()

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


### The estimation of the Gaussian parameters.
# INPUT:    binned histogram data
# OUTPUT:   estimated mean and standard deviation of the normal distribution

from scipy.stats import norm


# generate in-silico data
mean = 0.0
std  = 1.0
N    = 10000
X    = np.random.normal(mean, std, N)

density, bin_edges = np.histogram(X, density=True)
L_distr = np.diff(bin_edges) * density.cumsum()
R_distr = np.flip(np.diff(bin_edges), 0) * np.flip(density, 0).cumsum()

# I cannot use quantile of 1: it's infinity.
# nor can I use quantiles of 0: it's minus infinity
L_distr = L_distr[1:]
R_distr = R_distr[1:]
G = norm(loc=0, scale=1)
# I need corresponding quantile


z0 = bin_edges[1:-1]
z1 = np.flip(bin_edges, 0)[1:-1]

u0 = G.ppf(L_distr)
u1 = G.isf(R_distr)

u = np.concatenate((u0, u1))
z = np.concatenate((z0, z1))

u_sq  = sum(u**2)
u_sum = sum(u)
N     = len(L_distr) + len(R_distr)
z_sum = sum(z)
u_z   = sum(u * z)
det   = N * u_sq - u_sum**2

mean_lsq_estim = ( u_sq  * z_sum - u_sum * u_z) / det
sd_lsq_estim   = (-u_sum * z_sum + N     * u_z) / det
# this seriously sucks....



# check that the estimate of sigma ain't negative!!!





def least_squares_estimate():



# naive numpy weighted median



def numpy_weighted_median(data, weights=None):
    """Calculate the weighted median of an array/list using numpy."""
    import numpy as np
    if weights is None:
        return np.median(np.array(data).flatten())
    data, weights = np.array(data).flatten(), np.array(weights).flatten()
    if any(weights > 0):
        sorted_data, sorted_weights = map(np.array, zip(*sorted(zip(data, weights))))
        midpoint = 0.5 * sum(sorted_weights)
        if any(weights > midpoint):
            return (data[weights == np.max(weights)])[0]
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]
        if cumulative_weight[below_midpoint_index] - midpoint < sys.float_info.epsilon:
            return np.mean(sorted_data[below_midpoint_index:below_midpoint_index+2])
    return sorted_data[below_midpoint_index+1]





# implement the median estimation and standard deviation based on medi

# implement the estimation of the standard deviations.



# Wait, we don't want to use it directly.
# We might use it to spot too small and too big deviations, though.
# So we might invest more time to have something truly better.
# Like a spline.
# The small diffs are particularly troublesome.

# What we need to do really now, is to estimate the 
# standard deviations: the hereroscedasticity of the trend.
# Implement the estimator of the normal parameters based on 
# the historgram.

# Then, simply use it for passing to a function that will
# spread the IsoSpec infinitely-resolved peaks over the bins.

# Implement a binning based on one spectrum first.
# Then, think how to aggregate spectra.
# But this has to be done based on some similarity across runs.
# So it's more complicated: leave it for now.

# Instead, simply make a binning defined over what?
# The question is: should we have bins as large as
# the base peak groups, or smaller?
# If we assume, that more than one thing can explain them,
# then better smaller?
# It's a problem of continuous VS discrete.
