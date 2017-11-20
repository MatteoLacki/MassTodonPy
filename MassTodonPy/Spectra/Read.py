import numpy as np
from pyteomics import mzxml  # >= 3.41
from lxml import etree

from MassTodonPy.Spectra.Spectra import ExperimentalSpectrum as ExpSpec
from MassTodonPy.Spectra.Operations import round_n_trim


# TODO add assertions
def read_mzxml_spectrum(path,
                        mz_precision=10,
                        intensity_cut_off=0.0,
                        sum_intensity=False):
    """
    Generate a sequence of rounded and trimmed spectra from
    individual runs of the instrument.
    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    mz_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    intensity_cut_off : float
        The cut off value for peak intensity.
    sum_intensity : bool
        Report the total ion current.
    Returns
    -------
    out : generator
        Generates Spectra alone or together with their total ion counts.
    """
    with mzxml.read(path) as reader:
        for spectrum in reader:
            spectrum = ExpSpec(mz=spectrum['m/z array'],
                               intensity=spectrum['intensity array'])
            total_intensity = spectrum.intensity.sum()
            spectrum = ExpSpec(*round_n_trim(
                *spectrum,
                support_precision=mz_precision,
                value_cut_off=intensity_cut_off))
            if not sum_intensity:
                yield spectrum
            else:
                yield (spectrum, total_intensity)


def read_mzxml_spectrum_faster(path,
                               mz_precision=10,
                               intensity_cut_off=0.0,
                               sum_intensity=False):
    """
    Generate a sequence of rounded and trimmed spectra from
    individual runs of the instrument. A little bit faster.
    Parameters
    ----------
    path : str
        Path to the mzXml file containing the mass spectrum.
    mz_precision : integer
        The number of digits after which the support values get rounded.
        E.g. if set to 2, then number 3.141592 will be rounded to 3.14.
    intensity_cut_off : float
        The cut off value for peak intensity.
    sum_intensity : bool
        Report the total ion current.
    Returns
    -------
    out : generator
        Generates Spectra alone or together with their total ion counts.
    Remarks
    -------
        This code would not exist without the help of Joshua Klein.
        Thanks, Joshua!
    """
    for event, tag in etree.iterparse(path):
        if tag.tag.endswith("peaks"):
            spectrum = mzxml._decode_peaks(tag.attrib, tag.text)
            spectrum = ExpSpec(mz=spectrum['m/z array'],
                               intensity=spectrum['intensity array'])
            total_intensity = spectrum.intensity.sum()
            spectrum = ExpSpec(*round_n_trim(
                *spectrum,
                support_precision=mz_precision,
                value_cut_off=intensity_cut_off))
            if not sum_intensity:
                yield spectrum
            else:
                yield (spectrum, total_intensity)


def read_txt_spectrum(path,
                      mz_precision=10,
                      intensity_cut_off=0.0,
                      sum_intensity=False):
    """
    Read spectrum from a text file.
    Parameters
    ----------
    path : str
        Path to the *.txt text file containing the mass spectrum.
    sum_intensity : bool
        Report the total ion current.
    Returns
    -------
    out : Spectrum
    """
    mzs = []
    intensities = []
    with open(path) as f:
        for line in f:
            line = line.split()
            mzs.append(float(line[0]))
            intensities.append(float(line[1]))
    total_intensity = sum(intensities)
    spectrum = ExpSpec(*round_n_trim(np.array(mzs),
                                     np.array(intensities),
                                     mz_precision,
                                     intensity_cut_off))
    if not sum_intensity:
        return spectrum
    else:
        return (spectrum, total_intensity)
