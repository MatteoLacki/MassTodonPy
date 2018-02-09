from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet

from MassTodonPy.Misc.os import create_folder_if_needed


def plot_aggregated_precursors(masstodon,
                               path='aggregated_precusors.html',
                               mode='inline',
                               show_plot=True,
                               width=400,
                               height=400,
                               **kwds):
    """Plot intensity of precursors, neglecting the quenched charge.

    Parameters
    ==========
    intensities : list
        List of estimated intensities. Indices correspond to charge - 1.
    path : string
        Path to where to save the output html file.
        If not provided, MassTodon will use the folder it is called from.
    mode : string
        The mode of plotting a bokeh plot.
    width : integer
        The width of the plot.
    height : integer
        The height of the plot.

    """
    intensities = masstodon.report.get_aggregated_precursors()
    create_folder_if_needed(path)
    output_file(path, mode=mode)
    charges = list('q = {0}'.format(x) for x in range(1, 1+len(intensities)))
    p = figure(plot_width=400,
               plot_height=400,
               x_range=charges)
    p.xaxis.axis_label = 'Charge'
    p.yaxis.axis_label ='Estimated Intensity'
    data = ColumnDataSource({'intensity': intensities,
                             'charge': charges})
    bars = p.vbar(source=data, x='charge', top='intensity', width=1)
    hover_bars = HoverTool(renderers=[bars],
                           tooltips=[('charge', '@charge'),
                                     ('intensity', "@intensity{0,0}")],
                           mode='vline')
    p.add_tools(hover_bars)
    if show_plot:
        show(p)
    return p
