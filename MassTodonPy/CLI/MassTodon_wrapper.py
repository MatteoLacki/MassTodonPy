import json

from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_precursors import bokeh_aggregated_precursors
from MassTodonPy.Plot.bokeh_fragments_intensity import bokeh_fragments_intensity


def run_masstodon(masstodon_args,
                  output,
                  spectrum_plot_args={},
                  aggregated_precursors_plot_args={},
                  fragments_intensity_plot_args={},
                  _max_times_run_masstodon=10,
                  _verbose=False,
                  **kwds):
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
    while not Finished and i < _max_times_run_masstodon:
        try:
            masstodon = MassTodon(**masstodon_args)
            if _verbose:
                print('\tSaving results.')
            results_path = output + 'assigned_spectrum.csv'
            masstodon.report.write(results_path)
            masstodon.write(output)
            if _verbose:
                print('\tPlotting.')
            plot_path = output + 'assigned_spectrum.html'
            if spectrum_plot_args:
                spec = bokeh_spectrum(masstodon=masstodon,
                                      path=output+'assigned_spectrum.html',
                                      **spectrum_plot_args)
            if aggregated_precursors_plot_args:
                prec = bokeh_aggregated_precursors(masstodon=masstodon,
                                                   path=output+'aggregated_precusors.html',
                                                   **aggregated_precursors_plot_args)
            if fragments_intensity_plot_args:
                frag = bokeh_fragments_intensity(masstodon=masstodon,
                                                 path=output+'fragment_intensities.html',
                                                 **fragments_intensity_plot_args)
            Finished = True
        except ValueError:
            i += 1
            print('Failed '+str(i) +' times out of {0} already due to CVXOPT.'.format(_max_times_run_masstodon))
    if not Finished:
        print('We tried {0} times to run MassTodon, but CVXOPT always crushed.'.format(_max_times_run_masstodon))
        print('Switch to more stable Python2.7.')
