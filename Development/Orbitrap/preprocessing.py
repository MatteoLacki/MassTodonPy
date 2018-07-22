%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import numpy as np
import matplotlib.pyplot as plt
import os

from MassTodonPy.plotters.spectrum import plot_spectrum
from MassTodonPy.readers.from_npy import spectrum_from_npy

# mz_dummy        = np.array([1.1, 1.2, 1.3, 1.4, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0])
# intensity_dummy = np.array([0.1, 1.2, 2.3, 1.4, 0.2, 2.3, 4.4, 2.5, 1.6, 0.8, 2.9, 0.2])
# clusters_dummy = clusters(mz_dummy, intensity_dummy)
# plot(mz_dummy, intensity_dummy, list(clusters_dummy))


data_path = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

# cs = clusters(mz, intensity)
# colors_no = 10
# cmap = plt.get_cmap('tab10', colors_no)
# colors = [cmap(c % colors_no) for c in cs]
# plot(mz, intensity, colors)
# csa = clusters2array(mz, intensity)
# clusters_no = csa[-1] + 1


from MassTodonPy.Spectra.peak_clustering import mz_bitonic
import numpy as np


# around a second. That's something to rewrite later to C++.
clusters = mz_bitonic(mz, intensity, min_mz_diff=.15, abs_perc_dev=.2)




colors_no = 10
cmap = plt.get_cmap('tab10', colors_no)
colors = [cmap(c % colors_no) for c in clusters]
plot(mz, intensity, colors)
# it seems to work!!!
# final test: do we get a nice curve by looking at the diffs?

clusters = list(clusters)
mz_4_diff_plot = []
mz_diffs_4_diff_plot = []
for s, e in iter_clusters_from(clusters.__iter__()):
    mz_local = mz[s:e]
    mz_local_diffs = np.diff(mz_local)
    mz_4_diff_plot.extend(list(mz_local[0:-1]))
    mz_diffs_4_diff_plot.extend(list(mz_local_diffs))

plt.scatter(mz[:-1], np.diff(mz), c='grey')
plt.scatter(mz_4_diff_plot, mz_diffs_4_diff_plot, c='red')
plt.show()

from sklearn.linear_model import HuberRegressor

x = np.array(mz_4_diff_plot)
y = np.array(mz_diffs_4_diff_plot)
coefs = np.polyfit(x, y, 4)
parabola = np.poly1d(coefs)

x_parabola = np.linspace(min(x), max(x), 10000)
y_parabola = parabola(x_parabola)
plt.plot(x_parabola, y_parabola, color='red')
plt.scatter(x, y, s=1, c='grey')
plt.show()

# what is the median absolute error?
errors = y - parabola(x)
U = np.linspace(0,100, 100)/100.0
P = np.percentile(errors, 100*U)
plt.scatter(P, U)

plt.hist(parabola(x) - y)
plt.show()

# remove the observations that are too far away.
S = np.abs(errors) < np.percentile(np.abs(errors), 99.99)

plt.scatter(x, y, s=1, c=S)
plt.show()


new_parab = np.poly1d(np.polyfit(x[S], y[S], deg=4))
plt.plot(x_parabola, y_parabola, color='blue')
plt.plot(x_parabola, new_parab(x_parabola), color='red')
plt.scatter(x, y, s=1, c='grey')
plt.show()

# the impact of outliers seem truly negligable. Good.
# Leave order 4 polynomial, as is.
# Use it for other procedures.
# Forget about the robust fit: not so important.

# Order the code!
# The 4th order polynomial will be used to do what???

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
