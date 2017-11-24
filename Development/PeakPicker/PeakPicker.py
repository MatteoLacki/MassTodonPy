%load_ext autoreload
%autoreload 2

from math import ceil, log10
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.Precursor import Precursor
from MassTodonPy.MoleculeMaker.MoleculeMaker import get_molecules
from MassTodonPy.Spectra.Spectrum import Spectrum
from MassTodonPy.PeakPicker.PeakPicker import PeakPicker
from MassTodonPy.IsotopeCalculator.IsotopeCalculator import IsotopeCalculator

mol = get_dataset('substanceP')
# mol = get_dataset('ubiquitin')  # CORRUPTED!!!!

precursor_name = mol.precursor.name
precursor_fasta = mol.precursor.fasta
precursor_charge = mol.precursor.q
precursor_modifications = {10: {'C_carbo': {'H': 1,
                                            'N': 1,
                                            'O': -1}}}
blockedFragments = set([])
fragmentation_type = 'cz'
minimal_distance_between_charges = 5
m_over_z_precision = .05
spectrum_minimal_intensity = 100.
spectrum_percent_top_peaks = .95
_faster_mzxml = False

precursor = Precursor(name=precursor_name,
                      fasta=precursor_fasta,
                      q=precursor_charge,
                      modifications=precursor_modifications)

molecules = get_molecules([precursor],
                          blockedFragments,
                          fragmentation_type,
                          minimal_distance_between_charges)

spectrum = Spectrum(mol.spectrum,
                    3,
                    spectrum_minimal_intensity,
                    spectrum_percent_top_peaks,
                    _faster_mzxml)

list(spectrum)
spectrum.spectrum

isotope_calculator = IsotopeCalculator(mz_precision_digits=3,
                                       joint_probability=.999)

peak_picker = PeakPicker(molecules, isotope_calculator, mz_precision_digits)

spectrum.spectrum
peak_picker.represent_as_Graph(spectrum)
