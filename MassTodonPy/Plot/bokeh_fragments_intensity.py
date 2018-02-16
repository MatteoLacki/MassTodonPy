from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet
from bokeh.resources import CDN
from bokeh.embed import file_html

from MassTodonPy.Parsers.Paths import parse_path
from MassTodonPy.Misc.io import create_folder_if_needed

def bokeh_fragments_intensity(masstodon,
                              path="assigned_spectrum.html",
                              mode="inline",
                              show_plot=True,
                              width=0,
                              height=0,
                              **kwds):
    """Make a plot of the intensity of fragments.

    Parameters
    ==========
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
    create_folder_if_needed(path)
    # output_file(path, mode=mode)
    output_file(path, mode='cdn')

    afi = masstodon.report.aggregeted_fragment_intensities()
    afi['z_minus'] = [-v for v in afi['z']]

    if width and height:
        p = figure(plot_width=1000,
                   plot_height=400)
    else:
        p = figure()

    p.xaxis.axis_label = 'Cleavage Site'
    p.yaxis.axis_label = 'Estimated Intensity'
    afi['x'] = list(range(len(afi['z'])))
    data = ColumnDataSource(afi)
    bars_c = p.vbar(source=data, x='x', top='c', width=.8)
    hover_bars_c = HoverTool(renderers=[bars_c],
                             tooltips=[('name', '@c_name'),
                                       ('intensity', "@c{0,0}")])
    p.add_tools(hover_bars_c)
    bars_z = p.vbar(source=data, x='x', top='z_minus', width=.8, color='red')
    hover_bars_z = HoverTool(renderers=[bars_z],
                             tooltips=[('name', '@z_name'),
                                       ('intensity', "@z{0,0}")])
    p.add_tools(hover_bars_z)
    # labels_c = LabelSet(source=data, x='x', y='c')
    # p.add_layout(labels_c)
    # labels_z = LabelSet(source=data, x='x', y='z')
    # p.add_layout([labels_z, labels_c])
    if show_plot:
        show(p)
    else:
        with open(path, 'w') as f:
            f.write(file_html(p, CDN, path))
    return p
