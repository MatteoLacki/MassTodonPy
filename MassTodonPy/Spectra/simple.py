from   bisect            import bisect_left, bisect_right
from   collections       import  Counter
import matplotlib.pyplot as plt
from   math              import floor, log10
import numpy             as np

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
                                                         sd,\
                                                         skewness


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
                 mz              = np.array([]),
                 intensity       = np.array([]),
                 sort            = True,
                 drop_duplicates = True,
                 only_positive   = True,
                 bc              = None,
                 mdc             = None,
                 mz_diff_model   = None,
                 min_mz_diff_bc  = None,
                 min_mz_diff_mdc = None):
        """Initialize the Spectrum."""
        self._store_names = ('m/z', 'intensity')
        self.clusters     = None
        self.mz, self.intensity = dedup_sort(mz,
                                             intensity,
                                             drop_duplicates,
                                             sort)
        self.mz        = self.mz[self.intensity > 0]
        self.intensity = self.intensity[self.intensity > 0]
        # parameters for spectra spawns as subspectra
        self.bc  = bc
        self.mdc = mdc
        self.mz_diff_model   = mz_diff_model
        self.min_mz_diff_bc  = min_mz_diff_bc
        self.min_mz_diff_mdc = min_mz_diff_mdc

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

    @property
    def min_mz(self):
        return min(self.mz)

    @property
    def max_mz(self):
        return max(self.mz)

    @property
    def interval(self):
        return self.min_mz, self.max_mz
 
    def mean_mz(self):
        """Mean m/z weighted by intensities."""
        return mean(self.mz, self.intensity)

    def sd_mz(self):
        """Get standard deviation of m/z with probs induced by intensity."""
        return sd(self.mz, self.intensity)

    def skewness_mz(self):
        return skewness(self.mz, self.intensity)

    def total_intensity(self):
        return self.intensity.sum()

    def l1(self):
        """Get l1 norm."""
        return sum(self.intensity)

    def l2(self):
        """Get l2 norm."""
        return np.linalg.norm(self.intensity)

    def get_mz_digits(self):
        """Return the smallest m/z difference within the first bc cluster."""
        mz, _    = next(self.iter_bc_clusters())
        min_diff = np.diff(mz)[0]
        return abs(floor(log10(min_diff)))

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
        self.min_mz_diff_bc = min_mz_diff

    def min_mz_diff_clustering(self,
                               min_mz_diff = 1.1):
        if self.min_mz_diff_bc:
            assert min_mz_diff >= self.min_mz_diff_bc, "This can break up clusters: reduce min_mz_diff to levels lower than used for bitonic clustering."
        self.mdc = min_diff_clustering(x          = self.mz,
                                       min_x_diff = min_mz_diff)
        self.min_mz_diff_mdc = min_mz_diff

    def iter_clusters(self, clustering):
        """Iterate over consecutive cluster ends."""
        for s, e in iter_cluster_ends(clustering):
            yield self.mz[s:e], self.intensity[s:e]

    def iter_bc_clusters(self):
        return self.iter_clusters(self.bc)

    def iter_mdc_clusters(self):
        return self.iter_clusters(self.mdc)

    def iter_subspectra(self, clustering):
        for s, e in iter_cluster_ends(clustering):
            yield self.__class__(mz              = self.mz[s:e],
                                 intensity       = self.intensity[s:e],
                                 sort            = False,
                                 drop_duplicates = True,
                                 bc              = self.bc[s:e],
                                 mdc             = self.mdc[s:e],
                                 mz_diff_model   = self.mz_diff_model,
                                 min_mz_diff_bc  = self.min_mz_diff_bc,
                                 min_mz_diff_mdc = self.min_mz_diff_mdc)

    def iter_bc_subspectra(self):
        return self.iter_subspectra(self.bc)

    def iter_mdc_subspectra(self):
        return self.iter_subspectra(self.mdc)

    def mz_lefts_mz_diffs_in_clusters(self):
        """Get the left ends of histogramed data and lengths of bases of the bins."""
        # the total number of diffs within clusters
        clusters_no = self.bc[-1] - self.bc[0]
        diffs_no    = len(self.mz) - clusters_no - 1
        mz_lefts    = np.zeros(shape = (diffs_no,), dtype = float)
        mz_diffs    = mz_lefts.copy()
        i_ = _i = 0
        for s, e in iter_cluster_ends(self.bc):
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
                      knots_no      = 1000,
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
            x = np.linspace(self.min_mz(), self.max_mz(), knots_no)
            y = self.mz_diff_model(x)
            plt.plot(x, y, c='red')
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

    #TODO: this is a good method for any clustering.
    def get_bc_stats(self):
        means = []
        sds   = []
        skewnesses = []
        counts= []
        total_intensities = []
        mz_spreads = []
        for local_mz, local_intensity in self.iter_bc_clusters():
            mean_mz       = mean(local_mz, local_intensity)
            sd_mz         = sd(local_mz, local_intensity, mean_mz)
            skewnesses_mz = skewness(local_mz, local_intensity, mean_mz, sd_mz)
            means.append(mean_mz)
            sds.append(sd_mz)
            skewnesses.append(skewnesses_mz)
            counts.append(len(local_mz))
            total_intensities.append(sum(local_intensity))
            mz_spreads.append(max(local_mz) - min(local_mz))
        o = tuple(map(np.array, [means, sds, skewnesses, counts, total_intensities, mz_spreads]))
        return o

    def fit_sd_mz_model(self,
                        model = polynomial,
                        fit_to_most_frequent = True,
                       *model_args,
                      **model_kwds):
        means, sds, skewnesses, counts, total_intensities, spreads = self.get_bc_stats()
        if fit_to_most_frequent:
            cnts, freq       = list(zip(*Counter(counts).items()))
            self.sd_mz_c     = counts == cnts[np.argmax(freq)]
            self.sd_mz_model = model(means[self.sd_mz_c],
                                     sds[self.sd_mz_c])
        else:
            self.sd_mz_model = model(means, sds)

    def plot_sd_mz(self,
                   plt_style = 'dark_background',
                   show      = True):
        self.sd_mz_model.plot(plt_style     = plt_style,
                              scatter_color = self.sd_mz_c,
                              show          = show)

def spectrum(mz        = np.array([]),
             intensity = np.array([]),
             sort      = True):
    spec = Spectrum(mz, intensity, sort)
    return spec
