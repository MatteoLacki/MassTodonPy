import numpy as np
from collections import Counter

from MassTodonPy.models.polynomial import polynomial
from MassTodonPy.models.spline     import spline
from MassTodonPy.Spectra.orbitrap.peak_clustering import bitonic_clustering,\
                                                         iter_cluster_ends,\
                                                         min_diff_clustering
from MassTodonPy.stats.simple_normal_estimators   import mean,\
                                                         sd,\
                                                         skewness
class PeakClustering(object):
    def __init__(self):
        self.clusters = None

    def _iter_cluster_ends(self):
        yield from iter_cluster_ends(self.clusters)

    def __iter__(self):
        for s, e in self._iter_cluster_ends():
            yield self.x[s:e], self.w[s:e]

    def left_ends_and_diffs(self):
        """Get the left ends of histogramed data and lengths of bases of the bins."""
        # the total number of diffs within clusters
        clusters_no = self.clusters[-1] - self.clusters[0]
        diffs_no    = len(self.x) - clusters_no - 1
        lefts       = np.zeros(shape=(diffs_no,), dtype=float)
        diffs       = lefts.copy()
        i_ = _i = 0
        for s, e in self._iter_cluster_ends():
            x = self.x[s:e]
            diff = np.diff(x)
            _i += len(x) - 1
            lefts[i_:_i] = x[:-1]
            diffs[i_:_i] = diff
            i_ = _i
        return lefts, diffs



class MinDiff(PeakClustering):
    def fit(self, x, w, min_mz_diff = 1.1):
        self.x = x
        self.w = w
        self.clusters = min_diff_clustering(self.x,
                                            min_mz_diff)


def min_diff_clust(x, w, min_mz_diff = 1.1):
    mdc = MinDiff()
    mdc.fit(x, w, min_mz_diff)
    return mdc





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

    def stats(self, cut_one_point_intervals=True):
        min_mz = []
        max_mz = []
        means  = []
        sds    = []
        skewnesses = []
        counts     = []
        total_intensities = []
        mz_spreads        = []
        for local_mz, local_intensity in self:
            min_mz.append(min(local_mz))
            max_mz.append(max(local_mz))
            mean_mz       = mean(local_mz, local_intensity)
            sd_mz         = sd(local_mz, local_intensity, mean_mz)
            skewnesses_mz = skewness(local_mz, local_intensity, mean_mz, sd_mz)
            means.append(mean_mz)
            sds.append(sd_mz)
            skewnesses.append(skewnesses_mz)
            counts.append(len(local_mz))
            total_intensities.append(sum(local_intensity))
            mz_spreads.append(max(local_mz) - min(local_mz))
        o = tuple(map(np.array, [min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, mz_spreads]))
        if cut_one_point_intervals:
            OK = o[0] < o[1]
            o  = [x[OK] for x in o]
        return o

    def fit_diff_model(self, 
                       model=spline,
                      *model_args,
                     **model_kwds):
        lefts, diffs = self.left_ends_and_diffs()
        self.diff_model = model(lefts, diffs, *model_args, **model_kwds)

    def fit_sd_model(self,
                     model = polynomial,
                     fit_to_most_frequent = True,
                    *model_args,
                   **model_kwds):
        min_mz, max_mz, means, sds, skewnesses, counts, total_intensities, spreads = self.stats()
        if fit_to_most_frequent:
            cnts, freq    = list(zip(*Counter(counts).items()))
            self.sd_mz_c  = counts == cnts[np.argmax(freq)]
            self.sd_model = model(means[self.sd_mz_c],
                                  sds[self.sd_mz_c])
        else:
            self.sd_model = model(means, sds)

    # TODO: this should contain more class-specific code
    def plot_sd(self,
                plt_style = 'dark_background',
                show      = True):
        self.sd_model.plot(plt_style = plt_style,
                           show      = show)

def bitonic_clust(x, w, 
                  min_mz_diff=.15,
                  abs_perc_dev=.2,
                  model_diff=spline,
                  model_diff_args=[],
                  model_diff_kwds={},
                  fit_to_most_frequent = True,
                  model_sd=polynomial,
                  model_sd_args=[],
                  model_sd_kwds={}):
    bc = Bitonic()
    bc.fit(x, w, min_mz_diff, abs_perc_dev)
    bc.fit_diff_model(model_diff,
                     *model_diff_args,
                    **model_diff_kwds)
    bc.fit_sd_model(model_sd,
                    fit_to_most_frequent,
                   *model_sd_args,
                  **model_sd_kwds)
    return bc


