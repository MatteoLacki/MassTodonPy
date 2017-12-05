install.packages('rbokeh')
library(rbokeh)
library(tidyverse)


path_to_csv = '/Users/matteo/Documents/MassTodon/MassTodonPy/Development/Bokeh/H4-qb.csv'
D = read_table(path_to_csv, col_names = F)
colnames(D) = c('mz', 'intensity')
D



base_plot = 
  figure(width=900, height=500, title="Mass Spectrum") %>%
  ly_rect(
    data    = D,
    xleft   = mz - .05,
    xright  = mz + .05,
    ybottom = 0,
    ytop    = intensity,
    hover   = "m/z = @mz | intensity = @intensity",
    color   = 'black'
  ) %>% 
  x_axis(label = "m/z") %>%
  y_axis(label = "Intensity")
