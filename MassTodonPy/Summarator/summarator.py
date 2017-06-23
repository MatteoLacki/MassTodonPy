from collections import Counter

def summarize_results(  peakPicker_stats,
                        spectra,
                        raw_masstodon_res   ):
    '''Summarize the results of MassTodon.'''
    summary = Counter()
    used_E_total_intensity = peakPicker_stats['total intensity of experimental peaks paired with isotopologues']

    # the intensity of peaks outside any graph
    unused_E_total_intensity = spectra['trimmed intensity']+spectra['total intensity after trim'] - used_E_total_intensity

    trimmed_intensity = spectra['trimmed intensity']

    for r in raw_masstodon_res:
        summary['L1_error'] += r['L1_error']
        # summary['L2_error'] += r['L2_error']
        summary['underestimates']+= r['underestimates']
        summary['overestimates'] += r['overestimates']

        if r['status'] != 'optimal':
            summary['L1_error_nonoptimal'] += r['L1_error']
            # summary['L2_error_nonoptimal'] += r['L2_error']
            summary['underestimates_nonoptimal']+= r['underestimates']
            summary['overestimates_nonoptimal'] += r['overestimates']

    if spectra['original total intensity'] > 0.0:
        summary['L1_error/original_total_intensity']= (summary['L1_error']+unused_E_total_intensity)/spectra['original total intensity']
        # summary['L2_error/original_total_intensity']= (summary['L2_error']+unused_E_total_intensity)/spectra['original total intensity']

    if spectra['total intensity after trim'] > 0.0:
        summary['L1_error/total_intensity_after_trim'] = (summary['L1_error']+trimmed_intensity)/spectra['total intensity after trim']
        # summary['L2_error/total_intensity_after_trim'] = (summary['L2_error']+trimmed_intensity)/spectra['total intensity after trim']

    if used_E_total_intensity > 0.0:
        summary['L1_error_on_scooped_mz/used_E_total_intensity'] = summary['L1_error']/used_E_total_intensity
        summary['underestimates/total_intensity_after_trim'] = summary['underestimates']/used_E_total_intensity
        summary['overestimates/total_intensity_after_trim'] = summary['overestimates']/used_E_total_intensity

    return summary
