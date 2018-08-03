from    math    import  inf
import  numpy   as      np

from MassTodonPy.plotters.spectrum import plot_spectrum


# around a second. That's something to rewrite later to C++.
def bitonic_iterator(mz, intensity, min_mz_diff=.15):
    """Cluster based on bitonicity of the intensity.

    If three consecutive intensities form a dump,
    I_{i-1} > I_i < I_{i+1}, then 'i' is considered
    a new cluster.
    Also, if m_i/z_i - m_{i-1}/z_{i-1} > min_mz_diff,
    then 'i' is considered a new cluster.

    Parametes
    ---------
    mz: np.array
        The recorded m/z values.
    intensity: np.array
        The recorded intensities.
    min_mz_diff : float
        The minimal m/z difference that separates clusters.
    """
    # previous two intensities: set up guardians
    i__ = 1.0
    _i_ = 0.0
    c   = -1    # cluster no: the first cluster will get tag '0'
    m_  = -inf # previous m/z 
    for _m, __i in zip(mz, intensity):
        big_mz_diff = _m - m_ > min_mz_diff
        if (__i >= _i_ and _i_ < i__) or big_mz_diff:
            c += 1
        yield c
        if big_mz_diff:
            i__ = 0.0
        else:
            i__ = _i_
        _i_ = __i
        m_  = _m


def clusters2array(mz, intensity,
                   clustering=bitonic_iterator,
                   **clustering_kwds):
    """Write elements of cluster iterator to a numpy array.

    Parameters
    ----------
    mz: np.array
        The recorded m/z values.
    intensity: np.array
        The recorded intensities.
    clustering: iterator
        Iterator of clusters.
    *clustering args:
        Unnamed list of arguments for the iterator.
    **clustering_kwds:
        Dictionary of pairs argument-value for the iterator.

    Returns
    -------
        np.array : an array of cluster assignments.
    """
    return np.fromiter(clustering(mz, intensity, **clustering_kwds),
                       dtype=int, count=len(mz))


def list_of_clusters(mz, intensity,
                     clustering=bitonic_iterator,
                     *clustering_args,
                     **clustering_kwds):
    """Write elements of cluster iterator to a list.

    Parameters
    ----------
    mz: np.array
        The recorded m/z values.
    intensity: np.array
        The recorded intensities.
    clustering: iterator
        Iterator of clusters.
    *clustering args:
        Unnamed list of arguments for the iterator.
    **clustering_kwds:
        Dictionary of pairs argument-value for the iterator.

    Returns
    -------
        list : list of cluster assignments.
    """
    return list(clustering(mz, intensity, **clustering_kwds))


def iter_cluster_ends(assignments):
    """Iterate over left and right ends of subsequent clusters defined by assignments.

    Provides iteration over 'sparse' representation of the assignments to clusters,
    i.e. tuples composed of the index of the first element in cluster,
    index of the last element of cluster.
    Clusters are assumed to appear orderly one after another.

    Parameters
    ----------
        assignments : iter of int
            Iterator with clusters' assignents.
    Yields
    ------
        tuple: indices marking the beginning and the end of a cluster.
    """
    c_ = 0
    _  = next(assignments) # cluster: 0 for a fresh generator.
    i_ = 0
    _i = 1
    for _c in assignments: # next cluster
        if _c == c_ + 1:
            yield i_, _i # start and end of cluster c_
            c_ += 1
            i_ = _i
        _i += 1
    yield i_, _i # start and end of the final cluster


def fix_local_clustering(mz, s, e, c,
                         abs_perc_dev = .2,
                         verbose      = False):
    """Fix bitonic clustering.

    Parameters
    ----------
    mz :np.array
        The recorded m/z values.
    s : int
        The index of the start of the cluster of peaks.
    e :int
        The index of the end of a cluster of peaks.
    c :int
        The number of current cluster of peaks.
    abs_perc_dev : float
        How big m/z deviation is tolerable.
    verbose : logical
        That's for the devel mode, which is not in place anyway.

    Yields
    ------
    np.array: assignment into clusters.
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


def mz_bitonic(mz, intensity,
               min_mz_diff  = .15,
               abs_perc_dev = .2):
    """Cluster points based on consecutive m/z distances and the bitonicity of the intensities.

    Parameters
    ----------
    mz: np.array
        The recorded m/z values.
    intensity: np.array
        The recorded intensities.
    min_mz_diff : float
        The minimal m/z difference that separates clusters.

    Returns
    -------
    np.array of ints: assignments into consecutive clusters.
    """
    # an array to store the clusters
    clusters = np.full(shape = (len(mz),), fill_value = 0, dtype = int)
    # the clustering iterator
    bc = bitonic_iterator(mz, intensity, min_mz_diff = min_mz_diff)
    # cluster count
    c = 0
    for s, e in iter_cluster_ends(bc):
        lclust = fix_local_clustering(mz, s, e, c, abs_perc_dev)
        clusters[s:e] = lclust
        if lclust[0] != lclust[-1]:
            c = lclust[-1]
        else:
            c += 1
    return clusters



def show_clustering_is_good(data_path = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1',
                            mz_name   = 'mz.npy',
                            intensity_name = 'in.npy'):
    """Show if the bitonic clustering meets the standard.

    Parameters
    ----------
    data_path : str
        The path to the folder containing the files.
    mz_name : str
        The name of the file that contains m/z values.
    intensity_name : str
        The name of the file that contains the intensity values.

    """
    mz, intensity = spectrum_from_npy(data_path)
    clusters = mz_bitonic(mz, intensity, min_mz_diff=.15, abs_perc_dev=.2)
    idx = (1146.95 < mz) & (mz < 1149.1)
    plot_spectrum(mz[idx], intensity[idx], clusters[idx])
