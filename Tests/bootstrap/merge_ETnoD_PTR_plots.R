setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/')
source('get_data.R')
source('wave_height_150.R')
source('wave_velocity_300.R')

plot1 = cowplot::plot_grid(
  WH_150_anion_no_cation_plot_no_legend, 
  WH_150_estimates_plot,
  nrow=2,  align='v'
)

plot2 = cowplot::plot_grid(
  WV_300_anion_no_cation_plot_no_legend, 
  WV_300_estimates_plot,  
  nrow=2,  align='v'
)

plot3 = cowplot::plot_grid(
  legend, 
  plot1, 
  plot2,
  nrow = 3,
  rel_heights = c(.025,1,1)
)

date = format(Sys.time(), "%a_%b_%d_%X_%Y")

cowplot::ggsave(paste0('visual/ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''),plot3, limitsize = F, dpi=600, width = 210, height = 210, units = 'mm', scale=1)
cowplot::ggsave(paste0('visual/ETnoD_PTR_meetingCation___', date, '.pdf', sep='', collapse=''),plot3, limitsize = F, dpi=600, width = 210, height = 210, units = 'mm', scale=1)

