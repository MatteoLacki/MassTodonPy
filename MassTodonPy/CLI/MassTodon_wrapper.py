import json
import os
from pprint import pprint

from MassTodonPy.MassTodon import MassTodon
from MassTodonPy.CLI.PlotParser import plot_parser
from MassTodonPy.Plot.bokeh_spectrum import bokeh_spectrum
from MassTodonPy.Plot.bokeh_aggregated_fragments import bokeh_aggregated_fragments
from MassTodonPy.Plot.bokeh_aggregated_fragments_estimated import bokeh_aggregated_fragments_estimated
from MassTodonPy.Plot.bokeh_aggregated_precursors import bokeh_aggregated_precursors



def run_masstodon(args):
    '''Run MassTodonPy.

    Parameters
    ----------
    args : dict
        A dictionary with the parsed configuration of MassTodon and its plots.

    '''
    plots = plot_parser(args)      # Renaming variables for plotting
    if args['verbose']:
        print('Welcome to MassTodon!\n')
        print('Running MassTodon.')

    pprint(args)

    output = args.pop('output')
    if not os.path.exists(output):
        os.makedirs(output)

    max_times = args.pop("_max_times")

    i = 0
    Finished = False
    while not Finished and i < max_times:
        try:
            masstodon = MassTodon(**args)
            if args['verbose']:
                print('\tSaving results.')
            results_path = output + 'assigned_spectrum.csv'
            masstodon.report.write(results_path)
            masstodon.write(output)
            if args['verbose']:
                print('\tPlotting.')
            plot_path = output + 'assigned_spectrum.html'
            if plots['spectrum_plot_args']:
                spec = bokeh_spectrum(masstodon=masstodon,
                                      path=output+'assigned_spectrum.html',
                                      **plots['spectrum_plot_args'])
            if plots['aggregated_precursors_plot_args']:
                prec = bokeh_aggregated_precursors(masstodon=masstodon,
                                                   path=output+'aggregated_precusors.html',
                                                   **plots['aggregated_precursors_plot_args'])
            if plots['aggregated_fragments_plot_args']:
                frag = bokeh_aggregated_fragments(masstodon=masstodon,
                                                  path=output+'aggregated_fragments.html',
                                                  **plots['aggregated_fragments_plot_args'])
            if plots['estimated_aggregated_fragments_plot_args']:
                frag = bokeh_aggregated_fragments(masstodon=masstodon,
                                                  path=output+'estimated_aggregated_fragments.html',
                                                  **plots['estimated_aggregated_fragments_plot_args'])
            Finished = True
        except ValueError:
            i += 1
            print('Failed '+ str(i) +' times out of {0} already due to CVXOPT.'.format(_max_times_run_masstodon))

    if not Finished:
        print('We tried {0} times to run MassTodon, but CVXOPT always crushed.'.format(_max_times_run_masstodon))
        print('Switch to more stable Python2.7.')

    if args['verbose']:
        print('Thank you for using MassTodonPy!\n')
