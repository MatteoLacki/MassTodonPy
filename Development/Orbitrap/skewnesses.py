%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  Counter
import numpy                as      np
import matplotlib.pyplot    as      plt
from   time                 import  time
import intervaltree         as      iTree
import networkx             as      nx
from   networkx             import  connected_component_subgraphs
import pandas               as      pd

from MassTodonPy.readers.from_npy             import spectrum_from_npy
from MassTodonPy.Precursor.simple             import precursor
from MassTodonPy.IsotopeCalculator.simple     import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups import bitonic_clustering
from MassTodonPy.Spectra.simple               import spectrum
from MassTodonPy.Molecule.simple              import molecule
from MassTodonPy.Formula.Formula              import formula
from MassTodonPy.Molecule.simple              import Molecule
from MassTodonPy.stats.simple_normal_estimators   import mean,\
                                                         sd
from MassTodonPy.Data.Constants import infinity

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()
# TODO: replace this later on: a more complex model has to be fitted.
# and the estimation of the Standard Deviations of groups has to be 
# dependent upon more data points all across the spectrum.
spec.fit_sd_mz_model()
# spec.sd_mz_model.plot()

subspectra = list(spec.iter_mdc_subspectra())
spec

mz_local, _ = next(spec.iter_bc_clusters())
min_diff = np.diff(mz_local)[0]

# this shows, that it's worth exploring deviations from the symmetry around the average
# as potential measure of non-normality.
spec[1126.5:1126.8].plot()

means, sds, skewnesses, counts, total_intensities, mz_spreads = spec.get_bc_stats()
plt.style.use('dark_background')
plt.hist(skewnesses[skewnesses<infinity], bins=100)
plt.show()

finite_skew = skewnesses[skewnesses<infinity]
np.mean(finite_skew)
np.std(finite_skew)

good_skew = finite_skew[np.abs(finite_skew-np.median(finite_skew)) < 3 * np.std(finite_skew)]
plt.hist(good_skew)
plt.show()
len(skewnesses) - len(good_skew)
# this means that groups that are heavily skewed do not
# carry most of the probability.
plt.hist(total_intensities[skewnesses<infinity][np.abs(finite_skew-np.median(finite_skew)) < 3 * np.std(finite_skew)]/max(total_intensities))
plt.show()
# as with respect to the maximal intensity, theirs is fairly limited.

subspectra_bc = list(spec.iter_bc_subspectra())
subspectra_bc = np.array(subspectra_bc)
poor_groups = subspectra_bc[np.isin(counts, (13, 14))]
poor_groups[0].plot()
# again: deviation from center-symmetry
poor_groups[1].plot()

# DD = pd.DataFrame(dict(mean=means, sd=sds, intensity=total_intensities, cnt=counts, spread=mz_spreads))
# DD.to_csv(path_or_buf = '/Users/matteo/Projects/MassTodonPy/data/data_4_sd.csv', index=False)

