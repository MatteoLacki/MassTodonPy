%load_ext autoreload
%autoreload 2
%load_ext line_profiler

import numpy as np
import matplotlib.pyplot as plt
from time import time

from MassTodonPy.readers.from_npy    import spectrum_from_npy
from MassTodonPy.Precursor.Precursor import precursor
from MassTodonPy.plotters.spectrum import plot_spectrum
# from MassTodonPy.IsotopeCalculator.IsotopeCalculator import isotope_calculator
from MassTodonPy.IsotopeCalculator.simple import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups   import bitonic_clustering

data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)
bc            = bitonic_clustering(mz,
                                   intensity, 
                                   min_mz_diff  = .15,
                                   abs_perc_dev = .2)

fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name = "shit")
mols   = list(prec.molecules())


mol = mols[0]
mol.mean_mz
mol.monoisotopic_mz
f = mol.formula
list(f)

mol.formula

# Task: enhance the calculation
# t0 = time()
# isotopologues = [m.isotopologues() for m in mols]
# t1 = time()
# 89 secs.
# this is actually not so bad,
# but how much can we save on using mean and sd?
# OK, have to have it done in Numpy to have a speed-up

iso_calc = isotope_calculator()
iso_calc._masses
iso_calc._probabilities
iso_calc._mean_mass
iso_calc._mean_variance


iso_calc.monoisotopic(mol.formula)
iso_calc.mean(mol.formula, q=10, g=12)
iso_calc.sd(mol.formula,   q=10, g=12)

iso_calc(mol.formula)



from IsoSpecPy.IsoSpecPy import IsoSpec

isospec = IsoSpec((100,200),
                  ((1.1, 2.31), (12.4, 13.45)),
                  ((.9, .1), (.5,.5)),
                  .999)

%%timeit
isospec.getConfs()

IsoSpecify(str(mol.formula), .99)


masses_npy = np.array(list(masses))

def getConfsNumpyBuffered(isospec):
    ffi = FFI()
    m, l, c = isospec.getConfsRaw()
    obj_cnt = len(m)
    obj_size = np.dtype('float64').itemsize
    m = np.frombuffer(ffi.buffer(m, obj_cnt * obj_size),
                      dtype = np.dtype('float64'))
    l = np.frombuffer(ffi.buffer(l, obj_cnt * obj_size),
                      dtype = np.dtype('float64'))
    return m, l

def getConfsNumpy(isospec):
    m, l, c = isospec.getConfsRaw()
    return np.array(list(m)), np.array(list(l))

%%timeit
m, l = getConfsNumpy(isospec)

%%timeit
m, l = getConfsNumpyBuffered(isospec)

%%timeit
m, l, c = isospec.getConfsRaw()

%%timeit
m, l, c = isospec.getConfs()



%%timeit
getConfsNumpy(isospec)




%%timeit
means = [iso_calc.mean(m.formula) for m in mols]




iso_calc.get_envelope('C200H400')
iso_calc('C200H400')



atoms_in_all_formulas = set(prec.formula._storage.keys())



mol.formula
A = iso_calc(mol.formula)
masses, probs = A.atoms, A.masses


plot_spectrum(masses, probs, show=False)
plt.scatter(mol.mean_mz, 0, c= 'red')
plt.show()

