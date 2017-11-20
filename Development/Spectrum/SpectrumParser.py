%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Spectra.SpectrumParser import \
    read_n_preprocess_spectrum
from MassTodonPy.Spectra.operations import remove_lower_quantile


subP = get_dataset('substanceP')
spectrum_from_spectrum = read_n_preprocess_spectrum(subP.spectrum,
                                                    100., 1.0, 2)

spectrum_path = '/Users/matteo/Documents/MassTodon/Data/'
spectrum_txt = spectrum_path + 'subP_spectrum.txt'
spectrum_from_txt = read_n_preprocess_spectrum(spectrum_txt,
                                               100., 1.0, 2)

spectrum_mzxml = spectrum_path + 'FRL_220715_ubi_952_ETD_40ms_01.mzXML'
spectrum_from_mzxml = read_n_preprocess_spectrum(spectrum_mzxml,
                                                 100., 1.0, 2)

%%time
remove_lower_quantile(*spectrum_from_mzxml)

def remove_lower_quantile(mz, intensities, retained_percentage=.95):
    """
    Remove a portion of the smallest peaks that cover
    1-retained_percentage of the total ion current.
    Parameters
    ----------
    mz : array
        The m/z values.
    intensity : array
        The intensities to bin.
    retained_percentage : float
        The percentage of the original total ion current
        to be retained after trimming.
    Returns
    -------
    effective_cut_off : float
        The minimal peak height.
    """
    mz_res, intensities_res = tuple(
        np.array(x) for x in
        zip(*sorted(zip(mz, intensities), key=itemgetter(1))))

    i = 0
    S = 0.0
    total = intensities_res.sum()
    while True:
        S += intensities_res[i]/total
        if S < 1.0-retained_percentage:
            i += 1
        else:
            break
    effective_cut_off = intensities_res[i]
    return effective_cut_off
