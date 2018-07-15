%load_ext autoreload
%autoreload 2
%load_ext line_profiler



import numpy as np
import matplotlib.pyplot as plt
import os

mz          = np.array([1.1, 1.2, 1.3, 1.4, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 2.9, 3.0])
intensity   = np.array([0.1, 1.2, 2.3, 1.4, 0.2, 2.3, 4.4, 2.5, 1.6, 0.8, 2.9, 0.2])

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

bc = bitonic_clustering(mz, intensity)
clusters = list(bc)
plot(mz, intensity, clusters)


data_path = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = get_spectrum(data_path)
clusters = np.fromiter(bitonic_clustering(mz, intensity),
                       dtype=int,
                       count=len(mz))
colors_no = 10
cmap = plt.get_cmap('tab10', colors_no)
colors = [cmap(c % colors_no) for c in clusters]
plot(mz, intensity, colors)

# empirical test suite:
mz_test_idx = np.logical_and( mz >= 858.0, mz <= 859.4 )
mz_test = mz[mz_test_idx]
intensity_test = intensity[mz_test_idx]

plot(mz_test, intensity_test)

