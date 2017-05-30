library('ggplot2')
library('tidyverse')
library(jsonlite)

setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')

# list.files()
# i = 1
exp_set = read_json(path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/experimental_settings.json')
exp_set = exp_set %>% lapply(unlist) %>% sapply(function(x) x) %>% t


for( i in 0:51 ){
  setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')
  
  fitted_data = paste0(i,'_shortData.csv')
  origin_data = paste0(i,'_substanceP.csv')
  
  WH = exp_set[i+1,1]
  WV = exp_set[i+1,2]
  
  output_name = paste0(i,'_WH-',WH,'_WV-',WV,'.pdf')
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
  
  setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/spectra_plots/')
  
  ggsave(filename=output_name, plot=p, width=40, height=6)
  
}
