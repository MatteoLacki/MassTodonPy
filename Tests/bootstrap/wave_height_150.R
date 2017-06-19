# source('get_data.R')
get_algo = Vectorize(function(string) strsplit(string,'_')[[1]][1]) 
check_if_precursor = Vectorize(function(string) tail(strsplit(string,'_')[[1]],1) == 'precursor') 

WH_150 = DD %>% filter(WH==150)

ETnoD_cols = find_among_colnames('ETnoD', WH_150)

WH_150_ETnoD = WH_150[,c(info_columns, ETnoD_cols)]
WH_150_Errors= WH_150[,c(info_columns, error_columns)]
# WH_150_NoFrags = WH_150[,c(info_columns, nofrag_columns)]
WH_150_anion_no_cation = WH_150[,c(info_columns, anion_no_cation_columns)]
WH_150_anion_no_cation =  
  WH_150_anion_no_cation %>%
  gather_('algo', 'val', anion_no_cation_columns) %>%
  mutate(algo = get_algo(algo)) %>% split(.$status)

# WH_150_anion_no_reaction =  
#   WH_150[,c(info_columns, total_reaction_columns)] %>%
#   gather_('algo', 'val', anion_no_cation_columns) %>%
#   mutate(algo = get_algo(algo)) %>% split(.$status)



WH_150_ETnoD = 
  WH_150_ETnoD %>% 
  gather_('algo', 'val', ETnoD_cols) %>%
  filter(!check_if_precursor(algo))  %>%
  mutate(algo = get_algo(algo))

WH_150_ETnoD = WH_150_ETnoD %>% split(.$status)

WH_150_estimates_plot = 
  WH_150_ETnoD$bootstrap %>%
  ggplot(aes(x=factor(WV), y=val)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = WH_150_ETnoD$real, color = 'red', shape=19) +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1)
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  ) +
  scale_y_continuous(
    labels   = scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = 'PTR'
      )
  )+
  xlab('Wave Velocity (Wave Height set to 150)') + 
  ylab('ETnoD')


WH_150_anion_no_cation_plot = 
  WH_150_anion_no_cation$real %>%
  ggplot(aes(x=factor(WV), y=1-val, color=algo)) +
  theme_minimal() +
  geom_jitter(shape=19, height = 0.0, width = .25) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  scale_colour_discrete(
    name   = 'Pairing Algorithm: ',
    labels = c('Basic', 'Intermediate', 'Advanced')
  )+
  xlab('Wave Velocity (Wave Height set to 150)') +
  ylab('Reaction') + 
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

legend <- cowplot::get_legend(WH_150_anion_no_cation_plot)

WH_150_anion_no_cation_plot_no_legend = 
  WH_150_anion_no_cation$real %>%
  ggplot(aes(x=factor(WV), y=1-val, color=algo)) +
  theme_minimal() +
  geom_jitter(shape=19, height = 0.0, width = .25) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  xlab('Wave Velocity (Wave Height set to 150)') +
  ylab('Reaction') + 
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(color=FALSE)


cowplot::plot_grid(legend, WH_150_anion_no_cation_plot_no_legend, WH_150_estimates_plot, nrow=3, align='v', rel_heights = c(.1, 1,1))
