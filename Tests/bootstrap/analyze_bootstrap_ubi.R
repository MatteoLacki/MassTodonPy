library(tidyverse); library(stringr); library(car); library(jsonlite)

D = list.files('ubi_14_07_2017_csv', full.names = T) %>%
  lapply(read.csv) %>%
  lapply(as_data_frame) 

E = lapply(D, function(d) d %>% mutate(preActive = factor(preActive))) %>% bind_rows(.id='ID_2') 

# E %>% group_by( preActive, supActive, retentionTime, run, precursorMZ  ) %>%
#   summarize( n = n()  ) %>% data.frame()
# E %>% group_by( preActive, supActive, retentionTime, run ) %>%
#   summarize( n = n()  ) %>% filter( n > 251 )

initial_run_files = list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun') %>% tools::file_path_sans_ext()

E = 
  E %>% 
  mutate( files = tools::file_path_sans_ext(files) ) %>%
  split( ifelse(.$files %in% initial_run_files, 'initial_run', 'normal_run'), drop = T )

E$initial_run$run = 'initial_run'

DE = E$normal_run %>% mutate(Q = round(8500/precursorMZ)) %>% split(.$Q)



prep_data_4_plot = function(DataForPlot)
  DataForPlot %>%
  select( ID, retentionTime, preActive, supActive, real_or_bootstrap, intensity_within_tolerance, contains('prob.ETnoD'), contains('anion_approached_cation') ) %>%
  select( -contains('precursor') ) %>%
  select( -contains('precursor') ) %>%
  select( -contains('frag') ) %>%
  select( -intensity_within_tolerance ) %>%
  select( -contains('anion') ) %>%
  gather("prob_type", "prob_val", -c(1:5) ) %>%
  mutate(
    Algorithm = sapply( strsplit(prob_type, '[.]'), '[[', 1 )
  ) %>%
  filter( Algorithm == 'basic') %>%
  split(.$real_or_bootstrap)

make_nice_plot = function(data_4_plot)
  data_4_plot$boot %>%
  ggplot(aes(x=ordered(retentionTime), y=prob_val)) +
  theme_grey() +
  geom_boxplot() +
  geom_point(data = data_4_plot$real, color = 'red', shape=19) +
  facet_grid(supActive ~ preActive, labeller = label_both) +
  scale_y_continuous(
    limits = c(0,1),
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-.,
      labels= scales::percent,
      name  = "PTR"
    )
  ) +
  theme(legend.position="top")+
  theme_minimal()+
  xlab('Retention Time') +
  ylab('ETnoD')

estimates_plot_9 = DE$'9' %>% prep_data_4_plot %>% make_nice_plot()
estimates_plot_6 = DE$'6' %>% prep_data_4_plot %>% make_nice_plot()
estimates_plot = cowplot::plot_grid(estimates_plot_9, estimates_plot_6, nrow=2, align='v', labels = c('Q = 9', 'Q = 6'))
