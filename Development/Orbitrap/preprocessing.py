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
bc            = bitonic_clustering(mz, intensity, min_mz_diff=.15, abs_perc_dev=.2)
bc.fit_mz_diffs()
bc.plot_mz_diffs()
bc.fit_mz_diffs(model = polynomial)
bc.plot_mz_diffs()






# implement the finding of the erroneous cases.

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
