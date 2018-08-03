from bisect import bisect_left, bisect_right
from MassTodonPy.Measure.Measure    import Measure
from MassTodonPy.plotters.spectrum  import plot_spectrum



class OrbitrapSpectrum(Measure):
    def __init__(self,
                 mz,
                 intensity,
                 sort = True):
        self._store_names = ('m/z', 'intensity')
        self.mz = mz
        self.intensity = intensity
        if sort:
            self.sort()

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

    def plot(self, **kwds):
        plot_spectrum(self.mz, self.intensity, **kwds)

    def __getitem__(self, key):
        return self.mz[key], self.intensity[key]

    def zoom(self, mz_left, mz_right):
        id_s = bisect_left(  self.mz, mz_left)
        id_e = bisect_right( self.mz, mz_right)
        return self.mz[id_s:id_e], self.intensity[id_s:id_e]
