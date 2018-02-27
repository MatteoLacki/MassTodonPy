from bokeh.embed import file_html
from bokeh.models import HoverTool, LabelSet
from bokeh.plotting import ColumnDataSource, figure, output_file, show as show_plot
from bokeh.resources import CDN


from MassTodonPy.Misc.io import create_folder_if_needed
from MassTodonPy.Plot.Misc import aggregate_fragments


def bokeh_aggregated_fragments_estimated(masstodon,
                               path='aggregated_estimated_fragments.html',
                               mode='inline',
                               show=False,
                               width=0,
                               height=0,
                               _offset=.2,
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
    show : boolean
        Should we show the plot in your default browser?
    plot_width : integer
        The width of the plot.
    plot_height : integer
        The height of the plot.

    """
    data = aggregate_fragments(report=masstodon.report,
                               fasta=masstodon.precursor.fasta,
                               offset=_offset)
    source = ColumnDataSource(data=data)
    create_folder_if_needed(path)
    output_file(path, mode=mode)
    if not height or not width:
        p = figure()
    else:
        p = figure(plot_width=int(width),
                   plot_height=int(height))
    p.xaxis.visible = False
    p.yaxis.axis_label = "Intensity"
    intensity_bars = p.vbar(x='x', width=2*_offset, top='intensity',
                            color="red", source=source)
    aa_labels = LabelSet(x='x', y=0.0, text='aa',
                         level='glyph',
                         y_offset=-20,
                         # x_offset=-10,
                         source=source, render_mode='css')
    p.add_layout(aa_labels)
    hover_bars = HoverTool(renderers=[intensity_bars],
                           mode='vline',
                           tooltips=[('amino acid', '@aa'),
                                     ('intensity', "@intensity{0.}"),
                                     ('probability',"@probability{0.00}%")])
    p.add_tools(hover_bars)
    if show:
        show_plot(p)
    else:
        with open(path, 'w') as f:
            f.write(file_html(p, CDN, path))
    return p
