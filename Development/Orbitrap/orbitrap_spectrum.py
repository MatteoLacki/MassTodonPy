from MassTodonPy.Measure.Measure import Measure


class OrbitrapSpectrum(Measure):
    def __init__(self, mz, intensity, sort = True):
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

    def plot(self):
        

orbi_spec = OrbitrapSpectrum(mz, intensity)
orbi_spec.mz
orbi_spec.intensity