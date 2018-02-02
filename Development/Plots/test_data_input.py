from bokeh.plotting import ColumnDataSource, figure, output_file, show
from bokeh.models import HoverTool, Span, LabelSet

tools = "crosshair pan wheel_zoom box_zoom undo redo reset box_select save".split(" ")

for_fig = { 'plot_width': 100,
            'plot_height': 200,
            'tools': tools,
            'x_label': 'di[a]'}
p.line()
p = figure(**for_fig)
p.vbar(a=1)
