import json
from MassTodonPy.MassTodon import MassTodon

def perform_calculations(args, output):
    '''Run MassTodonPy.

    Parameters
    ----------
    args : dict
        A dictionary with the parsed configuration of MassTodon.
    output : str
        Path to the output.

    '''

    if args.verbose:
        print('Running MassTodon.')
    masstodon = MassTodon(**masstodon_args)

    if args.verbose:
        print('Saving results.')
    results_path = output + 'assigned_spectrum.csv'
    masstodon.report.write(results_path)

    if args.plot:
        if args.verbose:
            print('Plotting.')
        plot_path = output_path + 'assigned_spectrum.html'
        width = int(args.width) if args.width else None
        height = int(args.height) if args.height else None
        masstodon.report.plot(plot_path, width=width, height=height)
    print('\nThank you for using MassTodonPy!\nSupport us with your grant money!\n')
