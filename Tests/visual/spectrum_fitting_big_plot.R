library('ggplot2')
library('tidyverse')
library(jsonlite)

setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')

# list.files()
# i = 1
exp_set = read_json(path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/experimental_settings.json')
exp_set = exp_set %>% lapply(unlist) %>% sapply(function(x) x) %>% t

make_plots = function(i){  
  fitted_data = paste0(i,'_shortData.csv')
  origin_data = paste0(i,'_substanceP.csv')
  
  WH = exp_set[i+1,1]
  WV = exp_set[i+1,2]
  
  plot_title = paste0('WH = ',WH, ', WV = ',WV)
  
  D = read.csv(fitted_data) %>% tbl_df
  R = read.csv(origin_data) %>% tbl_df
  
  p = D %>% select(mz_L, mz_R, tot_estimate, tot_intensity) %>%
    rename( L = mz_L, R=mz_R ) %>%
    mutate( mean_mz = factor((L+R)/2) ) %>%
    gather( "tag", "value", 3:4 ) %>%
    ggplot()+
    geom_bar(aes(x=mean_mz, y=value, fill=tag), stat='identity', position = 'dodge') +
    theme_dark()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 100)+
    xlab('Average mz of the given tolerance interval')+
    ylab('Intensity')+
    labs(title =plot_title )
  p
}


plots = lapply(0:51, make_plots)

library(cowplot)

big_plot = plot_grid(plotlist=plots, nrow=52)
ggsave(filename='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/spectra_plots/big_plot.pdf', plot=big_plot, width=40, height=6*52, limitsize = F)




# D0 = read.csv(paste0(0,'_shortData.csv')) %>% tbl_df
# D2 = read.csv(paste0(2,'_shortData.csv')) %>% tbl_df
# D0
# D2
