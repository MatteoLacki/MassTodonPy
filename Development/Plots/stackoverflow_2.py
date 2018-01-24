from collections import defaultdict

import numpy as np
from scipy.stats import norm

from bokeh.plotting import show, figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.palettes import Viridis6

RT_x = np.linspace(118, 123, num=50)

mass_spec = defaultdict(list)
for scale, mz in [(1.0, 83), (0.9, 55), (0.6, 98), (0.4, 43), (0.2, 39), (0.12, 29)]:
    mass_spec["RT"].append(RT_x)
    mass_spec["RT_intensity"].append(norm(loc=120.4).pdf(RT_x) * scale)
    mass_spec['MZ_tip'].append(mz)
    mass_spec['Intensity_tip'].append(scale)
mass_spec['color'] = Viridis6

source = ColumnDataSource(mass_spec)

p = figure(plot_height=400)
p.multi_line(xs='RT', ys='RT_intensity', legend="Intensity_tip",
             line_width=5, line_color='color', line_alpha=0.6,
             hover_line_color='color', hover_line_alpha=1.0,
             source=source)

p.add_tools(HoverTool(show_arrow=False, line_policy='next', tooltips=[
    ('MZ', '@MZ_tip'),
    ('Rel Intensity', '@Intensity_tip')
]))

show(p)
