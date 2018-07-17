%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import numpy as np
import matplotlib.pyplot as plt
import os

mz_dummy        = np.array([1.1, 1.2, 1.3, 1.4, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0])
intensity_dummy = np.array([0.1, 1.2, 2.3, 1.4, 0.2, 2.3, 4.4, 2.5, 1.6, 0.8, 2.9, 0.2])

def sorted(x):
    return np.all(np.diff(mz_r))

def plot(mz, intensity, clusters=None, show=True):
    """Make a simple visualization of the data."""
    plt.vlines(x=mz, ymin=[0], ymax=intensity)
    if clusters:
        plt.scatter(x=mz, y=np.zeros(len(mz)), c=clusters)
    if show:
        plt.show()

def get_spectrum(data_path):
    return [ np.load(os.path.join(data_path, p)) for p in ('mz.npy', 'in.npy') ]

# task:
#   bitonic marker of modes
#   I don't remember what it's used for
#       ok, to estimate the whole lot of underlying gaussian parameters.

# the programme:
#   write and test (lineprof) pure python code for the task.
    # rewrite in C and expose using cffi
    # better even: rewrite in bloody Golang, becuase other things are there.
    # then we can even assume zeros existing.
    #   easier criterion: zero intensity triggers the end of old clustering
    # THERE ARE NO ZEROS THERE!!!


# min_mz_diff can be established statistically, as the first percentile of mass diffs.

def bitonic_clustering(mz, intensity, min_mz_diff=.15):
    """Generate cluster numbers based on bitonicity of the intensity."""
    # previous two intensities
    i__ = 1.0
    _i_ = 0.0
    # cluster no
    c = -1 # the first cluster will get tag '0'
    # previous m/z 
    m_ = - np.inf
    for _m, __i in zip(mz, intensity):
        big_mz_diff = _m - m_ > min_mz_diff
        if (__i - _i_ >= 0.0 and _i_ < i__) or big_mz_diff:
            c += 1
        yield c
        if big_mz_diff:
            i__ = 0.0
        else:
            i__ = _i_
        _i_ = __i
        m_  = _m

def get_clusters(clustering, count):
    return np.fromiter(clustering, dtype=int, count=count)

def clusters2array(mz, intensity, clustering=bitonic_clustering, **clustering_kwds):
    return np.fromiter(clustering(mz, intensity, **clustering_kwds),
                       dtype=int, count=len(mz))

def clusters(mz, intensity, clustering=bitonic_clustering, **clustering_kwds):
    return list(clustering(mz, intensity, **clustering_kwds))

clusters_dummy = clusters(mz_dummy, intensity_dummy)
plot(mz_dummy, intensity_dummy, list(clusters_dummy))

data_path = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = get_spectrum(data_path)
plot(mz, intensity)

cs = clusters(mz, intensity)

colors_no = 10
cmap = plt.get_cmap('tab10', colors_no)
colors = [cmap(c % colors_no) for c in cs]
plot(mz, intensity, colors)

csa = clusters2array(mz, intensity)
clusters_no = csa[-1] + 1


def iter_clusters_from(assignments):
    """ Get left and right ends of subsequent clusters defined by assignments.

    It iterates over 'sparse' representation of the assignments to clusters,
    i.e. tuples composed of the index of the first element in cluster,
    index of the last element of cluster.
    Clusters are assumed to appear orderly one after another.

    Args:
        assignments (iter of int): assignents to cluster.
    """
    c_ = next(assignments) # cluster: 0 for a fresh generator.
    i_ = 0
    _i = 1
    for _c in assignments: # next cluster
        if _c == c_ + 1:
            yield i_, _i # start and end of cluster c_
            c_ += 1
            i_ = _i
        _i += 1
    yield i_, _i # start and end of the final cluster


bc = bitonic_clustering(mz, intensity, min_mz_diff=.15)
s_e_clusters = iter_clusters_from(bc)


s, e = next(s_e_clusters)
mz_local, intensity_local = mz[s:e], intensity[s:e]


def get_local_clustering(mz, s, e, c,
                         abs_perc_dev = .2,
                         verbose = False):
    """Find out, if indices much be changed.

    Args:
        mz (np.array of floats): the recorded m/z values.
        s (int): the index of the start of the cluster of peaks.
        e (int): the index of the end of a cluster of peaks.
        c (int): the number of current cluster of peaks.
        abs_perc_dev (float):
    Yields:
        np.array of ints: assignment into clusters.
    """
    clusters = np.full(shape = (e-s,), fill_value = c, dtype=int)
    mz_local = mz[s:e]
    # differences of consecutive m/z values
    mz_diffs = np.diff(mz_local)
    # the median should be a stable values to compare to
    # as there are vastly more similar diffs than other.
    me_mz_diff = np.median(mz_diffs)
    cc = c
    signal_border = abs_perc_dev * me_mz_diff
    for i, mz_diff in enumerate(mz_diffs):
        if abs(mz_diff - me_mz_diff) > signal_border:
            cc += 1
            if verbose:
                print('Found poor clustering between {} and {}'.format(s,e))
        clusters[i+1] = cc
    return clusters


bc = bitonic_clustering(mz, intensity, min_mz_diff=.15)
s_e_clusters = iter_clusters_from(bc) # put this into loop definition
abs_perc_dev = .2
clusters = np.full(shape = (len(mz),),
                   fill_value = 0,
                   dtype = int)
c = 0
for s, e in s_e_clusters:
    lclust = get_local_clustering(mz, s, e, c, abs_perc_dev, True)
    clusters[s:e] = lclust
    if lclust[0] != lclust[-1]:
        c = lclust[-1]
    else:
        c += 1


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



# empirical test suite:
mz_test_idx = np.logical_and( mz >= 858.0, mz <= 859.4 )
mz_test = mz[mz_test_idx]
intensity_test = intensity[mz_test_idx]

# plot(mz_test, intensity_test, )
# zero_intensity = intensity_test == 0.0
# plt.scatter(mz_test[zero_intensity], intensity_test[zero_intensity])
# plt.show()

# investigating clustering based on m/z distances
cs = clusters(mz_test, intensity_test, min_mz_diff=.6)
colors = [cmap(c % colors_no) for c in cs]
plot(mz_test, intensity_test, colors)


# global robust regression approach:
    # there should be ~ 6 to 8 times more small diffs than big ones


x = mz[0:-1]
y = np.diff(mz)
plt.scatter(x**2, y)
plt.show()

plt.hist(np.log(y) - 2 * np.log(x), bins=100, density=True)
plt.show()

n, bins, patches = plt.hist(x, 50, density=True, facecolor='g', alpha=0.75)

# good idea: fix the local assignment!


