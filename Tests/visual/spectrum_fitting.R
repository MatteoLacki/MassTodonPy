library('ggplot2')
library('tidyverse')
library(jsonlite)

setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')

# list.files()
# i = 0
exp_set = read_json(path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/experimental_settings.json')
exp_set = exp_set %>% lapply(unlist) %>% sapply(function(x) x) %>% t

for( i in 0:51 ){
  setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/')
  
  fitted_data = paste0(i,'_shortData.csv')
  
  WH = exp_set[i+1,1]
  WV = exp_set[i+1,2]
  
  output_name = paste0(i,'_WH-',WH,'_WV-',WV,'.pdf')
  plot_title = paste0('WH = ',WH, ', WV = ',WV)
  
  D = read.csv(fitted_data) %>% tbl_df %>% 
      select(mz_L, mz_R, tot_estimate, tot_intensity, where) %>%
      rename( L = mz_L, R=mz_R ) %>%
      mutate( mean_mz = factor((L+R)/2) )
  
  D = D %>% gather( "tag", "value", 3:4 ) %>%
    mutate( tag = ifelse(where == 'not_explainable', 'atheoretic', ifelse(tag=='tot_estimate', 'Estimated', 'Observed') ) )
  
  p = D %>%
    ggplot()+
    geom_bar(aes(x=mean_mz, y=value, fill=tag), stat='identity', position = 'dodge') +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 100)+
    xlab('Average mz of the given tolerance interval')+
    ylab('Intensity')+
    labs(title =plot_title )+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  setwd('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/spectra_plots/')
  
  info_size = nrow(D)
    
  ggsave(filename=output_name, plot=p, width=100*info_size/890, height=6, limitsize = F)
}
