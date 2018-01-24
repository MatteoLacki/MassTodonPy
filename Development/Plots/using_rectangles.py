from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.palettes import viridis
# from bokeh.themes.theme import Theme
#
# theme = Theme(filename="./theme.yaml")

output_file('pure_rectangles.html')
TOOLS = "crosshair pan wheel_zoom box_zoom undo redo reset box_select lasso_select save".split(" ")

p = figure(plot_width=400, plot_height=400, tools=TOOLS)

data = {'charge':   [20, 10, 30, 32],
        'color':    viridis(4),
        'name':     ['A', 'B', 'C', 'D'],
        'top':      [1, 2, 3, 4],
        'bottom':   [0, 1, 2, 3],
        'left':     [0, 1, 2, 3],
        'right':    [1, 2, 3, 4],
        'L':        [.4, 1.4, 2.4, 3.4],
        'R':        [.6, 1.6, 2.6, 3.6]}

data2= {'name':     ['E'],
        'top':      [5],
        'bottom':   [4],
        'left':     [4],
        'right':    [5],
        'L':        [4.4],
        'R':        [4.6]}

source = ColumnDataSource(data)
source2= ColumnDataSource(data2)

fat_rectangles = p.quad(top='top', bottom='bottom',
                        left='left', right='right',
                        color='color',
                        source=source)
slim_rectangles = p.quad(top='top', bottom='bottom', left='L', right='R',
                    source=source, color='black')

slim_rectangles2= p.quad(top='top', bottom='bottom', left='L', right='R',
                    source=source2, color='black')

hover = HoverTool(renderers=[fat_rectangles],
                  tooltips=[('Yani', "@name"),
                            ('Matteo', "@charge")])
p.add_tools(hover)
show(p)
