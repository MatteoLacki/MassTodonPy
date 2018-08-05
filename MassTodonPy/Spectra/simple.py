from   bisect import bisect_left, bisect_right
import matplotlib.pyplot as plt
import numpy  as np

from MassTodonPy.arrays.operations  import dedup_sort
from MassTodonPy.Data.Constants     import eps, infinity
from MassTodonPy.Measure.Measure    import Measure
from MassTodonPy.models.spline      import spline
from MassTodonPy.models.polynomial  import polynomial
from MassTodonPy.plotters.spectrum  import plot_spectrum


from MassTodonPy.Spectra.orbitrap.peak_clustering import bitonic_clustering,\
                                                         iter_cluster_ends,\
                                                         min_diff_clustering

from MassTodonPy.stats.simple_normal_estimators   import mean,\
                                                         sd


class Spectrum(Measure):
    """Prepare experimental spectrum.

    Parameters
    ----------
    spectrum : string, tuple, list, or Spectrum
        The path can end up with two extension:
        * txt, for spectra saved in tab separated format.
        * mzXml, for spectra saved with mxXml format.
        The tuple or list consist of two numpy arrays,
        one with m over z ratios and the other with intensities.
    mz_digits : float or int
        The number of digits after which the floats get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    min_intensity : float
        Experimental peaks with lower height will be trimmed.
    percent_top_peaks : float
        Percentage of the heighest peaks in the spectrum to be included.
    """
    def __init__(self,
                 mz        = np.array([]),
                 intensity = np.array([]),
                 sort      = True,
                 drop_duplicates = True,
                 bc              = None,
                 mdc             = None,
                 mz_diff_model   = None):
        """Initialize the Spectrum."""
        self._store_names = ('m/z', 'intensity')
        self.clusters     = None
        self.mz, self.intensity = dedup_sort(mz,
                                             intensity,
                                             drop_duplicates,
                                             sort)
        # parameters for spectra spawns as subspectra
        if bc is not None:
            self.bc = bc
        if mdc is not None:
            self.mdc = mdc
        if mz_diff_model is not None:
            self.mz_diff_model = mz_diff_model

    @property
    def mz(self):
        """Get mass over charge ratios"""
        return self.atoms

    @mz.setter
    def mz(self, mz):
        """Set m/z ratios."""
        self.atoms = mz

    @property
    def intensity(self):
        """Get intensities."""
        return self.masses

    @intensity.setter
    def intensity(self, intensity):
        """Set intensities."""
        self.masses = intensity

    def __getitem__(self, interval):
        id_s = bisect_left(self.mz, interval.start)
        id_e = bisect_right(self.mz, interval.stop)
        return self.__class__(self.mz[id_s:id_e], 
                              self.intensity[id_s:id_e],
                              sort = False)

    def mean_mz(self):
        """Mean m/z weighted by intensities."""
        return mean(self.mz, self.intensity)

    def sd_mz(self):
        """Get standard deviation of m/z with probs induced by intensity."""
        return sd(self.mz, self.intensity)

    def total_intensity(self):
        return self.intensity.sum()

    def l1(self):
        """Get l1 norm."""
        return sum(self.intensity)

    def l2(self):
        """Get l2 norm."""
        return np.linalg.norm(self.intensity)

    def trim_intensity(self, cut_off):
        """Trim intensities below the provided cut off.

        Parameters
        ----------
        cut_off : float

        """
        self.trim(cut_off)

    def plot(self, 
             plt_style = 'dark_background',
             show      = True,
             clusters  = 'bc'):
        clusters = self.bc if clusters == 'bc' else self.mdc
        plot_spectrum(mz        = self.mz,
                      intensity = self.intensity,
                      clusters  = clusters,
                      plt_style = plt_style,
                      show      = show)

    def bitonic_clustering(self,
                           min_mz_diff  = .15,
                           abs_perc_dev = .2):
        self.bc = bitonic_clustering(x = self.mz,
                                     w = self.intensity,
                                     min_x_diff   = min_mz_diff,
                                     abs_perc_dev = abs_perc_dev)

    def min_mz_diff_clustering(self,
                               min_mz_diff = 1.1):
        self.mdc = min_diff_clustering(x          = self.mz,
                                       min_x_diff = min_mz_diff)

    def iter_clusters(self, clustering):
        """Iterate over consecutive cluster ends."""
        for s, e in iter_cluster_ends(np.nditer(clustering)):
            yield self.mz[s:e], self.intensity[s:e]

    def iter_bc_clusters(self):
        return self.iter_clusters(self.bc)

    def iter_mdc_clusters(self):
        return self.iter_clusters(self.mdc)

    def iter_subspectra(self, clustering):
        for s, e in iter_cluster_ends(np.nditer(clustering)):
            yield self.__class__(mz        = self.mz[s:e],
                                 intensity = self.intensity[s:e],
                                 sort      = False,
                                 drop_duplicates = True,
                                 bc        = self.bc[s:e],
                                 mdc       = self.mdc[s:e],
                                 mz_diff_model = self.mz_diff_model)

    def iter_bc_subspectra(self):
        return self.iter_subspectra(self.bc)

    def iter_mdc_subspectra(self):
        return self.iter_subspectra(self.mdc)

    def mz_lefts_mz_diffs_in_clusters(self):
        """Get the left ends of histogramed data and lengths of bases of the bins."""
        # the total number of diffs within clusters
        diffs_no = len(self.mz) - self.bc[-1] - 1
        mz_lefts = np.zeros(shape = (diffs_no,), dtype = float)
        mz_diffs = mz_lefts.copy()
        i_ = _i = 0
        for s, e in iter_cluster_ends(np.nditer(self.bc)):
            mz      = self.mz[s:e]
            mz_diff = np.diff(mz)
            _i += len(mz) - 1
            mz_lefts[i_:_i] = mz[:-1]
            mz_diffs[i_:_i] = mz_diff
            i_ = _i
        return mz_lefts, mz_diffs

    def fit_mz_diff_model(self, 
                          model=spline,
                         *model_args,
                        **model_kwds):
        mz_lefts, mz_diffs = self.mz_lefts_mz_diffs_in_clusters()
        self.mz_diff_model = model(mz_lefts, mz_diffs, *model_args, **model_kwds)


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
            mz_lefts, mz_diffs = self.mz_lefts_mz_diffs_in_clusters()
            try:
                signal = self.mz_diff_model.is_signal
                c = np.array(['papayawhip' if s else 'yellow' for s in signal])
            except AttributeError:
                c = 'papayawhip'
            plt.scatter(mz_lefts,
                        mz_diffs,
                        c = c,
                        s = 1.5)
        if show and (all_diffs or cluster_diffs):
            plt.show()



def spectrum(mz        = np.array([]),
             intensity = np.array([]),
             sort      = True):
    spec = Spectrum(mz, intensity, sort)
    return spec
