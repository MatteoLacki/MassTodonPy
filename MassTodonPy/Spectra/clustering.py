from MassTodonPy.Spectra.orbitrap.peak_clustering import bitonic_clustering,\
                                                         iter_cluster_ends,\
                                                         min_diff_clustering
class PeakClustering(object):
    def __init__(self):
        self.clusters = None

    def _iter_cluster_ends(self):
        yield from iter_cluster_ends(self.clusters)

    def __iter__(self):
        for s, e in self._iter_cluster_ends():
            yield self.x[s:e], self.w[s:e]


class MinMZdiff(PeakClustering):
    def fit(self, x, w, min_mz_diff = 1.1):
        self.x = x
        self.w = w
        self.clusters = min_diff_clustering(self.x,
                                            min_mz_diff)


class Bitonic(PeakClustering):
    def __iter__(self):
        for s, e in self._iter_cluster_ends():
            yield self.x[s:e], self.w[s:e]

    def fit(self, x, w, 
                  min_mz_diff  = .15,
                  abs_perc_dev = .2):
        self.x = x
        self.w = w
        self.clusters = bitonic_clustering(self.x,
                                           self.w,
                                           min_mz_diff,
                                           abs_perc_dev)

