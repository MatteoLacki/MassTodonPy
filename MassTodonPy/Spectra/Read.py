import numpy as np
from pyteomics import mzml  # >= 3.41
from pyteomics import mzxml  # >= 3.41

from MassTodonPy.Data.Constants import infinity
from MassTodonPy.Data.Constants import eps
from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Spectra.ExperimentalSpectrum import ExperimentalSpectrum

# TODO add checks about the MS number
def read_mz_file(path,
                 mz_precision=10,
                 intensity_cut_off=eps,
                 format='mzxml'):
    """Read mzXML spectra.

    Generate a sequence of rounded and trimmed spectra from
    individual runs of the instrument.
    Parameters
    ----------
    path : str
        Path to the mass spectrum file (*.mzxml, *.mzml, *.txt).
    mz_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    intensity_cut_off : float
        The cut off value for peak intensity.
    format : str
        Either 'mzxml' or 'mzml'.
    Returns
    -------
    out : ExperimentalSpectrum
    """
    format = format.lower()
    reader = {'mzxml': mzxml, 'mzml': mzml}[format]
    with reader.read(path) as info:
        for spectrum in info:
            spectrum = ExperimentalSpectrum(spectrum['m/z array'],
                                            spectrum['intensity array'])
            spectrum.trim_intensity(intensity_cut_off)
            spectrum.round_mz(mz_precision)
            yield spectrum


def read_txt_file(path,
                  mz_precision=10,
                  intensity_cut_off=eps):
    """Read spectrum from a text file.

    Parameters
    ----------
    path : str
        Path to the mass spectrum file (*.mzxml, *.mzml, *.txt).
    mz_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    intensity_cut_off : float
        The cut off value for peak intensity.
    Returns
    -------
    out : ExperimentalSpectrum
    """
    mzs = []
    intensities = []
    with open(path) as f:
        for line in f:
            line = line.split()
            mzs.append(float(line[0]))
            intensities.append(float(line[1]))
    spectrum = ExperimentalSpectrum(mzs, intensities)
    spectrum.trim_intensity(intensity_cut_off)
    spectrum.round_mz(mz_precision)
    return spectrum


def read_spectrum(path='',
                  mz_precision=10,
                  intensity_cut_off=eps):
    """Read mzXML spectra.

    Parameters
    ----------
    path : str
        Path to the mass spectrum file (*.mzxml, *.mzml, *.txt).
    mz_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    intensity_cut_off : float
        The cut off value for peak intensity.
    Returns
    -------
    out : ExperimentalSpectrum

    """
    assert path is not '', "Provide a spectrum to analyze!"
    file_path, file_name, file_ext = parse_path(path)
    file_ext = file_ext.lower()
    if file_ext in ('.txt', ''):
        spectrum = read_txt_file(path,
                                 mz_precision,
                                 intensity_cut_off)
    elif file_ext in ('.mzml', '.mzxml'):
        spectra = read_mz_file(path,
                               mz_precision,
                               intensity_cut_off,
                               file_ext[1:])
        spectrum = sum(spectra)
    else:
        raise NotImplementedError
    return spectrum
