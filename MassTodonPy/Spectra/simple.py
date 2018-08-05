from   bisect import bisect_left, bisect_right
import numpy  as     np

from MassTodonPy.arrays.operations  import dedup_sort
from MassTodonPy.Data.Constants     import eps, infinity
from MassTodonPy.Measure.Measure    import Measure
from MassTodonPy.plotters.spectrum  import plot_spectrum


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
                 drop_duplicates = True):
        """Initialize the Spectrum."""
        self._store_names = ('m/z', 'intensity')
        self.mz, self.intensity = dedup_sort(mz,
                                             intensity,
                                             drop_duplicates,
                                             sort)

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
             show      = True):
        plot_spectrum(mz        = self.mz,
                      intensity = self.intensity,
                      clusters  = self.clusters(),
                      plt_style = plt_style,
                      show      = show)

    def cluster(self, clustering_algorithm):
        self.clustering = clustering_algorithm(self.mz, self.intensity)


def spectrum(mz        = np.array([]),
             intensity = np.array([]),
             sort      = True):
    spec = Spectrum(mz, intensity, sort)
    return spec