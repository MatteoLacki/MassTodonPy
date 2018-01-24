from bokeh.core.properties import value
from bokeh.io import show, output_file
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
from bokeh.models import HoverTool

output_file("stacked.html", mode='inline')

fruits = [1,2, 5,6,7, 10]
vegies = [21, 23, 25, 26, 29, 39]
years = ["2015", "2016", "2017"]
years2 = ["201", "206", "207"]
colors = ["#c9d9d3", "#718dbf", "#e84d60"]
colors2 = ["#718dbf", "#e84d60", "#c9d9d3"]
widths = [.1, .2, .3, .4, .5, .6]

data1 = {'fruits' : fruits,
        'widths' : widths,
        '2015'   : [2, 1, 4, 0, 2, 4],
        '2016'   : [5, 3, 4, 2, 4, 6],
        '2017'   : [3, 2, 4, 4, 5, 3]}

data2 = {'vegies' : vegies,
        'widths' : widths,
        '2015'   : [2, 1, 4, 0, 2, 4],
        '2016'   : [5, 3, 4, 2, 4, 6],
        '2017'   : [3, 2, 4, 4, 5, 3]}

source1 = ColumnDataSource(data=data1)
source2 = ColumnDataSource(data=data2)

TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")
p = figure(plot_height=250, title="Fruit Counts by Year", tools=TOOLS)

fruits = p.vbar_stack(stackers=years, x='fruits', width='widths',
                      color=colors, source=source1)

# vegies = p.vbar_stack(stackers=years, x='vegies', width='widths',
#                       color=colors2, source=source2)

help(HoverTool)

# hover_invisible = HoverTool()
hover_invisible = HoverTool(tooltips=[('x', "@fruits{0,0}"), ("@height")])

# hover_invisible = HoverTool(renderers=[fruits],
#                             tooltips=[('x', "@x{0,0}"),
#                                       ('width', "@width{0,0}")],
#                             mode='vline')

p.add_tools(hover_invisible)
p.y_range.start = 0
p.x_range.range_padding = 0.1
p.xgrid.grid_line_color = None
p.axis.minor_tick_line_color = None
p.outline_line_color = None

show(p)
