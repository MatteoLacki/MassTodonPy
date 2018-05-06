import numpy as np
from numpy.random import multinomial, normal
from MassTodonPy.Spectra.operations import merge_runs


def makeRandomSpectrum(isotope_calc,
                       formulas,
                       quantities,
                       massDeviation,
                       jointProbability):
    """Simulate a mixture of isotopic envelopes.

    Parameters
    ----------
    isotopeCalc: IsotopeCalculator
        An isotope calculator instance.
    formulas : list of named.tuples
        e.g. Formula(name='precursor', formula='C63H98N18O13S1', q=3, g=0)
    quantities : list
        A list of total intensities of each molecular species to simulate from.
    massDeviation : float
        The standard deviation of the masses of isotopologues -
        theoretical equivalent of the mass resolution.
    jointProbability : float
        The joint probability of the theoretical isotopic envelope.

    Returns
    -------
    spectrum : tuple
        A tuple containing the theoretical spectrum:
        mass over charge values and intensities.
    """
    x0 = sum(quantities)

    def get_intensity_measure(formulas, quantities):
        for mol, quant in zip(formulas, quantities):
            ave_mz, ave_intensity = isotope_calc.isoEnvelope(
                atomCnt_str=mol.formula,
                jointProbability=jointProbability,
                q=mol.q, g=mol.g, prec_digits=2)
            ave_intensity = quant * ave_intensity
            yield ave_mz, ave_intensity

    mz_average, intensity = reduce(merge_runs,
                                   get_intensity_measure(mols,
                                                         quants))
    probs = intensity/sum(intensity)
    counts = np.array(multinomial(x0, probs), dtype='int')

    if sigma>0.0:
        spectrum = Counter()
        for m_average,cnt in zip(mz_average, counts):
            if cnt > 0:
                m_over_z = np.round(normal(loc=m_average, scale=sigma, size=cnt), prec_digits)
                spectrum.update(m_over_z)

        spectrum = np.array(spectrum.keys()), np.array([ float(spectrum[k]) for k in spectrum ])
    else:
        spectrum = (mz_average, counts)
    return spectrum
