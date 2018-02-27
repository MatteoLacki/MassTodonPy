def plot_parser(args):
    """Plot the arguments for plots."""
    plots = dict(spectrum_plot_args={},
                 aggregated_precursors_plot_args={},
                 aggregated_fragments_plot_args={},
                 estimated_aggregated_fragments_plot_args={})
    if args.pop('spectrum_plot'):
        plots['spectrum_plot_args'] = dict(width=args.pop('width_spectrum_plot'),
                                           height=args.pop('height_spectrum_plot'),
                                           show=args.pop('show_spectrum_plot'))
    if args.pop('aggregated_precursors_plot'):
        plots['aggregated_precursors_plot_args'] = \
            dict(width=args.pop('width_aggregated_precursors_plot'),
                 height=args.pop('height_aggregated_precursors_plot'),
                 show=args.pop('show_aggregated_precursors_plot'))
    if args.pop('aggregated_fragments_plot'):
        plots['aggregated_fragments_plot_args'] = \
            dict(width=args.pop('width_aggregated_fragments_plot'),
                 height=args.pop('height_aggregated_fragments_plot'),
                 show=args.pop('show_aggregated_fragments_plot'))
    if args.pop('estimated_aggregated_fragments_plot'):
        plots['estimated_aggregated_fragments_plot_args'] = \
            dict(width=args.pop('width_estimated_aggrageted_fragments_plot'),
                 height=args.pop('height_estimated_aggregated_fragments_plot'),
                 show=args.pop('show_estimated_aggregated_fragments_plot'))
    return plots
