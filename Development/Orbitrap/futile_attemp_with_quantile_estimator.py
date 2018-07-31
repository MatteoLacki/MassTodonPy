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
