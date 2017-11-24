from math import ceil, log10
from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.MoleculeMaker.Precursor import Precursor
from MassTodonPy.MoleculeMaker.MoleculeMaker import get_molecules
from MassTodonPy.Spectra.Spectrum import Spectrum


mol = get_dataset('substanceP')
mol = get_dataset('ubiquitin')

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
spectrum = mol.spectrum
spectrum_percent_top_peaks = .95
_faster_mzxml = False

mz_precision_digits = int(ceil(-log10(m_over_z_precision)))+1

precursor = Precursor(name=precursor_name,
                      fasta=precursor_fasta,
                      q=precursor_charge,
                      modifications=precursor_modifications)

molecules = get_molecules([precursor],
                          blockedFragments,
                          fragmentation_type,
                          minimal_distance_between_charges)

spectrum = Spectrum(spectrum,
                    mz_precision_digits,
                    spectrum_minimal_intensity,
                    spectrum_percent_top_peaks,
                    _faster_mzxml)

spectrum.total_intensity
spectrum.intensity_after_trim
