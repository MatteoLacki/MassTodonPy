import json
from MassTodonPy.MassTodon import MassTodon

def run_masstodon(args,
                  output_path,
                  _max_times=10,
                  verbose=False):
    '''Run MassTodonPy.

    Parameters
    ----------
    args : dict
        A dictionary with the parsed configuration of MassTodon and its plots.
    output_path : str
        Path to the output.
    verbose : boolean
        Call MassTodonPy in a verbose mode.

    '''
    i = 0
    Finished = False
    while i < _max_times or Finished:
        try:
            masstodon = MassTodon(**args)
            if args.verbose:
                print('Saving results.')
            results_path = output_path + 'assigned_spectrum.csv'
            masstodon.report.write(results_path)
            masstodon.write(output_path)
            if verbose:
                print('Plotting.')
            plot_path = output_path + 'assigned_spectrum.html'
            if args.spectrum_plot:
                spec = bokeh_spectrum(masstodon=masstodon,
                                      path=output_path + 'assigned_spectrum.html',
                                      **args)
            if aggregated_precusors_plot_args:
                prec = bokeh_aggregated_precursors(masstodon=masstodon,
                                                   path=output_path + 'aggregated_precusors.html',
                                                   **aggregated_precusors_plot_args)
            if fragments_intensity_plot_args:
                frag = bokeh_fragments_intensity(masstodon=masstodon,
                                                 path=output_path + 'fragment_intensities.html',
                                                 **fragments_intensity_plot_args)
            print('\nThank you for using MassTodonPy!\n')
            Finished = True
        except ValueError:
            i += 1
            print('Failed '+str(i) +' times out of {0} already due to CVXOPT.'.format(_max_times))
    if not Finished:
        print('We tried {0} times to run MassTodon, but CVXOPT always crushed.'.format(_max_times))
        print('Switch to more stable Python2.7.')
