# source('get_data.R')

get_algo = Vectorize(function(string) strsplit(string,'_')[[1]][1]) 
check_if_precursor = Vectorize(function(string) tail(strsplit(string,'_')[[1]],1) == 'precursor') 



WV_300 = DD %>% filter(WV==300)

ETnoD_cols = find_among_colnames('ETnoD', WV_300)

WV_300_ETnoD = WV_300[,c(info_columns, ETnoD_cols)]
WV_300_Errors= WV_300[,c(info_columns, error_columns)]

WV_300_anion_no_cation = WV_300[,c(info_columns, anion_no_cation_columns)]
WV_300_anion_no_cation =  
  WV_300_anion_no_cation %>%
  gather_('algo', 'val', anion_no_cation_columns) %>%
  mutate(algo = get_algo(algo)) %>% split(.$status)

WV_300_ETnoD = 
  WV_300_ETnoD %>% 
  gather_('algo', 'val', ETnoD_cols) %>%
  filter(!check_if_precursor(algo))  %>%
  mutate(algo = get_algo(algo))

WV_300_ETnoD = WV_300_ETnoD %>% split(.$status)

WV_300_estimates_plot = 
  WV_300_ETnoD$bootstrap %>%
  ggplot(aes(x=factor(WH), y=val)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = WV_300_ETnoD$real, color = 'red', shape=19) +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(
    labels   = scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = 'PTR'
    )
  )+
  xlab('Wave Height (Wave Velocity set to 300)') + 
  ylab('ETnoD')


WV_300_anion_no_cation_plot = 
  WV_300_anion_no_cation$real %>%
  ggplot(aes(x=factor(WH), y=1-val, color=algo)) +
  theme_minimal() +
  geom_jitter(shape=19, height = 0.0, width = .1) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  xlab('Wave Height (Wave Velocity set to 300)') +
  ylab('Reaction') + 
  theme(legend.position="top")
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

WV_300_anion_no_cation_plot_no_legend = 
  WV_300_anion_no_cation$real %>%
  ggplot(aes(x=factor(WH), y=1-val, color=algo)) +
  theme_minimal() +
  geom_jitter(shape=19, height = 0.0, width = .1) + 
  scale_x_discrete(position = "top") +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  xlab('Wave Height (Wave Velocity set to 300)') +
  ylab('Reaction') + 
  theme(legend.position="top")+
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1)
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  ) +
  guides(color=FALSE)


cowplot::plot_grid(WV_300_anion_no_cation_plot, WV_300_estimates_plot, nrow=2, align = 'v')
