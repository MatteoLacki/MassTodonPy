%load_ext autoreload
%autoreload 2
%load_ext line_profiler

from collections            import  Counter
import numpy                as      np
import matplotlib.pyplot    as      plt
from   time                 import  time
import intervaltree         as      iTree
import networkx             as      nx
from networkx               import  connected_component_subgraphs

from MassTodonPy.readers.from_npy             import spectrum_from_npy
from MassTodonPy.Precursor.simple             import precursor
from MassTodonPy.IsotopeCalculator.simple     import isotope_calculator
from MassTodonPy.Spectra.orbitrap.peak_groups import bitonic_clustering
from MassTodonPy.Spectra.simple               import spectrum
from MassTodonPy.Molecule.simple              import molecule
from MassTodonPy.Formula.Formula              import formula
from MassTodonPy.Molecule.simple              import Molecule

# generating subspectra
data_path     = '/Users/matteo/Projects/review_masstodon/data/PXD001845/numpy_files/20141202_AMB_pBora_PLK_10x_40MeOH_1FA_OT_120k_10uscans_928_ETciD_8ms_15SA_19precZ/1'
mz, intensity = spectrum_from_npy(data_path)

spec = spectrum(mz, intensity)
spec.bitonic_clustering()
spec.fit_mz_diff_model()
spec.min_mz_diff_clustering()
subspectra = list(spec.iter_mdc_subspectra())

# generating formulas
fasta  = "GAASMMGDVKESKMQITPETPGRIPVLNPFESPSDYSNLHEQTLASPSVFKSTKLPTPGKFRWSIDQLAVINPVEIDPEDIHRQALYLSHSRIDKDVEDKRQKAIEEFFTKDVIVPSPWTDHEGKQLSQCHSSKCTNINSDSPVGKKLTIHSEKSD"
charge = 24
prec   = precursor(fasta, charge, name = "")
mols   = list(prec.molecules())

# # meeting ends :)
# ends = []
# bc_groups = []
# for subspec in subspectra:
#     for bc_group in subspec.iter_bc_subspectra():
#         s, e = bc_group.interval
#         if s < e:
#             ends.append(s)
#             ends.append(e)
#             bc_groups.append(bc_group)
# # 13.5 sec
# # we are trying to figure out, if intervals are overlapping.
# mol_intervals = np.zeros(dtype=int,      shape=(2*len(mols),) )
# long_mols     = np.empty(dtype=Molecule, shape=(2*len(mols),) )
# for i, mol in enumerate(mols):
#     mol_intervals[2*i:2*(i+1)] = mol.interval(std_cnt=3)
#     long_mols[2*i:2*(i+1)]     = mol
# # 2.6
# sorting       = np.argsort(mol_intervals)
# mol_intervals = mol_intervals[sorting]
# long_mols     = long_mols[sorting]



# # division into subproblems.
# subproblems = []

# def bc_groups(subspectra)
#     for subspec in subspectra:
#         for bc_group in subspec.iter_bc_subspectra():
#             yield bc_group


# divmod(9, 2)

# encoding:
# -2 : start of empirical group
# -1 : end of empirical group
# divmod(K, 2): div: number of theoretical group
              # mod:
                # start = 0, end = 1



mol_intervals = np.zeros(dtype=float, shape=(2*len(mols),))
mol_indices   = np.zeros(dtype=int,   shape=(2*len(mols),))
for i, mol in enumerate(mols):
    mol_intervals[2*i:2*(i+1)] = mol.interval(std_cnt=3)
    mol_indices[2*i:2*(i+1)]   = i
is_mols_1 = np.full(shape      = mol_indices.shape,
                    fill_value = True)

i = 0
ends = []
bc_groups = []
subspec_indices = []
for subspec in subspectra:
    for bc_group in subspec.iter_bc_subspectra():
        s, e = bc_group.interval
        ends.append(s)
        ends.append(e)
        bc_groups.append(bc_group)
        subspec_indices.append(i)
        subspec_indices.append(i)
        i += 1

# these can be later reused for the in-depth validation
ends            = np.array(ends)
bc_groups       = np.array(bc_groups)
subspec_indices = np.array(subspec_indices)
is_mols_2       = np.full(shape      = subspec_indices.shape,
                          fill_value = False)

# merging into common lists.
all_indices = np.concatenate((mol_indices, subspec_indices))
is_mols     = np.concatenate((is_mols_1, is_mols_2))
all_mz      = np.concatenate((mol_intervals, ends))
sorting     = np.argsort(all_mz)
all_indices = list(all_indices[sorting])
is_mols     = list(is_mols[sorting])
all_mz      = list(all_mz[sorting])

# now I iterate over all m/z values of both signals and noise.
subproblems    = [] # a list of tuples of sets: mols and theory
open_theoretic = set([])
open_empiric   = set([])


for mz, is_mol, idx in zip(all_mz, is_mols, all_indices):
    t = T[i] # type of point.




subproblems.append(problem)
