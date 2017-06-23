from collections import Counter

def summarize_results(  spectra,
                        raw_masstodon_res   ):
    '''Summarize the results of MassTodon.'''
    summary = Counter()
    used_E_total_intensity = spectra['intensity of peaks paired with isotopologues']

    # the intensity of peaks outside any graph
    unused_E_total_intensity = spectra['trimmed intensity']+spectra['total intensity after trim'] - used_E_total_intensity
    trimmed_intensity = spectra['trimmed intensity']
    summary['intensity_original'] = spectra['original total intensity']
    summary['intensity_after_trim'] = spectra['total intensity after trim']
    summary['intensity_trimed_out'] = spectra['original total intensity'] - spectra['total intensity after trim']
    summary['intensity_within_tolerance'] = used_E_total_intensity

    for r in raw_masstodon_res:
        summary['L1_error'] += r['L1_error']
        summary['underestimates']+= r['underestimates']
        summary['overestimates'] += r['overestimates']

        if r['status'] != 'optimal':
            summary['L1_error_nonoptimal'] += r['L1_error']
            summary['underestimates_nonoptimal']+= r['underestimates']
            summary['overestimates_nonoptimal'] += r['overestimates']

    if spectra['original total intensity'] > 0.0:
        summary['L1_fit_error_and_unused_intensity/original_total_intensity']= (summary['L1_error']+unused_E_total_intensity)/spectra['original total intensity']

    if used_E_total_intensity > 0.0:
        summary['L1_error/intensity_within_tolerance'] = summary['L1_error']/used_E_total_intensity
        summary['underestimates/intensity_within_tolerance'] = summary['underestimates']/used_E_total_intensity
        summary['overestimates/intensity_within_tolerance'] = summary['overestimates']/used_E_total_intensity

    return summary
