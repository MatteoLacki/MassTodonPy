"""TODO: document properly."""

import matplotlib.pyplot as plt
import numpy as np

from MassTodonPy.arrays.operations          import dedup_sort
from MassTodonPy.models.spline              import spline
from MassTodonPy.plotters.spectrum          import plot_spectrum
from MassTodonPy.Spectra.peak_clustering    import mz_bitonic
from MassTodonPy.Spectra.peak_clustering    import iter_cluster_ends


class Cluster(object):
    def __init__(self, *_cluster_args, **_cluster_kwds):
        """Pass the arguments for the clustering algorithm here.

        In future, this will ease up copy contructors and similar things.
        The solution is taken from SKLearn.
        """
        raise NotImplementedError


    def cluster(self, mz, intensity, 
                drop_duplicates = True,
                sort            = True):
        self.mz, self.intensity = dedup_sort(mz, intensity, drop_duplicates, sort)
        self._cluster()


    def _cluster(self):
        raise NotImplementedError("This should be setting 'self.clusters'.")


    def __len__(self):
        """Get the number of clusters.

        Returns
        -------
            out : int
            The number of clusters.
        """
        return self.clusters[-1] + 1


    def iter_cluster_ends(self):
        """Iterate over consecutive cluster ends."""
        for i_, _i in iter_cluster_ends(np.nditer(self.clusters)):
            yield i_, _i


    def iter_mz_and_mz_diffs(self):
        """Iterate over m/z and m/z diffs.

        Iterate over tuples of arrays consisting of
        the left ends of spaces between consecutive
        peaks in a sequence of clusters.

        Yields:
        mz, mz_diff: tuples of np.arrays with mz values of all but the largest m/z in a cluster and the values of spaces between consecutive peaks.
        """
        for s, e in self.iter_cluster_ends():
            mz      = self.mz[s:e]
            mz_diff = np.diff(mz)
            yield mz[:-1], mz_diff


    def mz_and_mz_diffs(self):
        """Get m/z and m/z diffs."""
        diffs_no = len(self.mz) - len(self)
        mz = np.full(shape      = (diffs_no,),
                     fill_value = 0.0,
                     dtype      = float)
        mz_diffs = mz.copy()
        i_ = _i = 0
        for mz_l, mz_diff_l in self.iter_mz_and_mz_diffs():
            _i += len(mz_l)
            mz[i_:_i]       = mz_l
            mz_diffs[i_:_i] = mz_diff_l
            i_ = _i
        return mz, mz_diffs


    def fit_mz_diffs(self, model=spline, 
                     *model_args, **model_kwds):
        """Fit a spline to (m/z, Δm/z)."""
        mz_lefts, mz_diffs = self.mz_and_mz_diffs()
        self.mz_diff_model = model(mz_lefts, mz_diffs, *model_args, **model_kwds)


    def mz_diff(self, mz):
        """Return the fitted mz_diff for a given values of mz.

        Parameters
        ----------
        mz : np.array
            m/z values for which you need the values of the estimated m/z differences.
        """
        return self.mz_diff_model(mz)


    def plot(self,
             plt_style = 'dark_background',
             colors_no = 10,
             show      = True):
        """Plot the clustering on the mass spectrum.

        Parameters
        ----------
            plt_style : str
                The type of the matplotlib style used in the plot.
            colors_no : int
                The number of colors to cycle through.
            show : logical
                Immediately show the plot? Alternatively, just add it to the present canvas.

        """
        plot_spectrum(self.mz,
                      self.intensity,
                      self.clusters,
                      plt_style,
                      colors_no,
                      show)


    def plot_mz_diffs(self,
                      plt_style     = 'dark_background',
                      all_diffs     = True,
                      cluster_diffs = True,
                      trend         = True,
                      show          = True):
        """Plot the clustering on the mass spectrum.

        Parameters
        ----------
            plt_style : str
                The type of the matplotlib style used in the plot.
            all_diffs : logical
                Plot all the differences, or only the ones in clusters?
            cluster_diffs : logical
                Plot the differences that belong to the cluster.
            trend : logical
                Plot the fitted trendline.
            show : logical
                Immediately show the plot? Alternatively, just add it to the present canvas.

        """
        plt.style.use(plt_style)
        if all_diffs:
            # plot all the subsequent m/z differences as function of m/z
            plt.scatter(self.mz[:-1],
                        np.diff(self.mz),
                        c = 'blue',
                        s = .5)
        if trend:
            self.mz_diff_model.plot(plot_data = False,
                                    show      = False)
        if cluster_diffs:
            # mask those that are within clusters by bigger red dots.
            mz_lefts, mz_diffs = self.mz_and_mz_diffs()
            plt.scatter(mz_lefts,
                        mz_diffs,
                        c = 'papayawhip',
                        s = 1.5)
        if show and (all_diffs or cluster_diffs):
            plt.show()



class BitonicCluster(Cluster):
    """Clustering based on bitonic intensities."""
    def __init__(self, min_mz_diff = .15, abs_perc_dev = .2):
        self.min_mz_diff  = min_mz_diff
        self.abs_perc_dev = abs_perc_dev

    def _cluster(self):
        self.clusters = mz_bitonic(self.mz,
                                   self.intensity,
                                   self.min_mz_diff,
                                   self.abs_perc_dev)


def bitonic_clustering(mz, intensity,
                       min_mz_diff  = .15,
                       abs_perc_dev = .20):
    bc = BitonicCluster(min_mz_diff, abs_perc_dev)
    bc.cluster(mz, intensity, drop_duplicates = True, sort = True)
    return bc