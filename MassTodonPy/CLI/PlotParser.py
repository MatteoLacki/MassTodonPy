def plot_parser(args):
    """Plot the arguments for plots."""
    plots = dict(spectrum_plot_args={},
                 aggregated_precursors_plot={},
                 fragments_intensity_plot_args={})
    if args.pop('spectrum_plot'):
        plots['spectrum_plot_args'] = dict(width=args.pop('width_spectrum_plot'),
                                           height=args.pop('height_spectrum_plot'),
                                           show=args.pop('show_spectrum_plot'))
    if args.pop('aggregated_precursors_plot'):
        plots['aggregated_precursors_plot_args'] = dict(width=args.pop('width_aggregated_precursors_plot'),
                                                        height=args.pop('height_aggregated_precursors_plot'),
                                                        show=args.pop('show_aggregated_precursors_plot'))
    if args.pop('fragments_intensity_plot'):
        plots['fragments_intensity_plot_args'] = dict(width=args.pop('width_fragments_intensity_plot'),
                                                      height=args.pop('height_fragments_intensity_plot'),
                                                      show=args.pop('show_fragments_intensity_plot'))
    return plots
