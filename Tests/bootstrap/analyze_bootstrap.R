library('tidyverse'); library('stringr'); library('car')


D = list.files('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/CSV_03_07_2017_mzPrec-065', full.names = T) %>%
  lapply(read.csv) %>% 
  lapply(as_data_frame) %>% 
  bind_rows() %>%
  mutate(WH = factor(WH), WV = factor(WV) )


D150 = D %>% filter(WH == '150') %>% 
  select( contains('prob.ETnoD'), contains('anion_approached_cation'), real_or_bootstrap, WH, WV, ID, intensity_within_tolerance ) %>%
  gather("prob_type", "prob_val", 1:6 ) %>%
  split(.$real_or_bootstrap)


estimates_plot150 = 
  D150$boot %>%
  ggplot(aes(x=WV, y=prob_val)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = D150$real, color = 'red', shape=19) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1)
    # axis.title.x=element_blank()
    # axis.text.x=element_blank()
  ) +
  scale_y_continuous(
    labels   = scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = 'PTR'
    )
  )+
  xlab('Wave Velocity (Wave Height set to 1,5)') + 
  ylab('ETnoD')

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
get_algo = Vectorize(function(x) capFirst(strsplit(x,'[.]')[[1]][1]))

reaction_D150 = D150$boot %>%
  gather('algo','prob_reaction',1:3) %>%
  mutate(Algorithm = get_algo(algo)) 

reaction_plot150 =
  reaction_D150 %>%    
  ggplot(aes(x=WV, y=prob_reaction, fill=Algorithm, color=Algorithm)) + 
  geom_boxplot()+
  theme_minimal() +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  scale_x_discrete(position = "top") +
  xlab('Wave Velocity (Wave Height set to 1,5)') +
  ylab('Reaction') + 
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  
WH150 = cowplot::plot_grid(reaction_plot, estimates_plot, nrow=2, align='v')

################################################################################################

D300 = D %>% filter(WV == '300') %>% 
  select( contains('prob.ETnoD'), contains('anion_approached_cation'), real_or_bootstrap, WH, WV, ID, intensity_within_tolerance ) %>%
  gather("prob_type", "prob_val", 1:6 ) %>%
  split(.$real_or_bootstrap)


estimates_plot300 = 
  D300$boot %>%
  ggplot(aes(x=WH, y=prob_val)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = D300$real, color = 'red', shape=19) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1)
         # axis.title.x=element_blank()
         # axis.text.x=element_blank()
  ) +
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

reaction_D300 = D300$boot %>%
  gather('algo','prob_reaction',1:3) %>%
  mutate(Algorithm = get_algo(algo)) 

reaction_plot300 =
  reaction_D300 %>%    
  ggplot(aes(x=WH, y=prob_reaction, fill=Algorithm, color=Algorithm)) + 
  geom_boxplot()+
  theme_minimal() +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-., 
      labels=scales::percent,
      name  = "No Reaction"
    )
  )+
  scale_x_discrete(position = "top") +
  xlab('Wave Height (Wave Velocity set to 300)') +
  ylab('Reaction') + 
  theme(legend.position="top")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

WV300 = cowplot::plot_grid(reaction_plot300, estimates_plot300, nrow=2, align='v')


################################################################################################

reaction_prob_VS_etnod_precision_plot = 
  D %>% 
  select( contains('prob.ETnoD'), contains('anion_approached_cation'), real_or_bootstrap, WH, WV, ID, intensity_within_tolerance ) %>%
  gather("prob_type", "prob_val", 1:6 ) %>%
  gather("algo", "prob_reaction", 1:3)  %>%
  drop_na %>% group_by(ID, prob_type, WH) %>% 
  mutate( sd_estim = sd(prob_val), 
          var_estim = var(prob_val), 
          ave_intensity_within_tol = mean(intensity_within_tolerance), 
          ave_anion_approached_cation = mean(prob_reaction)) %>%
  ungroup %>% select(-real_or_bootstrap, -prob_type, -WH, -WV) %>%
  mutate(Algorithm = get_algo(algo)) %>%
  ggplot(aes(x=prob_reaction, y=1.0/var_estim, color=Algorithm, fill=Algorithm)) + 
  geom_point() +
  # geom_hex(bins=100) +
  geom_smooth(method='lm', formula=y ~ poly(x, 2)+poly(1/(1-x),4)) +
  theme_minimal() + 
  scale_x_continuous(
    labels=scales::percent,
    breaks=seq(0,1,.05)
  )+
  xlab('Average Probability of a Reaction') + 
  ylab('Precision of ETnoD Estimates')+
  theme(legend.position="bottom")


################################################################################################

WH_WV = cowplot::plot_grid(WH150, WV300, ncol=2, align='h', rel_widths = c(1.5,.75))

publication_plot = cowplot::plot_grid(WH_WV, reaction_prob_VS_etnod_precision_plot, nrow=2, rel_heights = c(1,.5))
# publication_plot = cowplot::plot_grid(WH_WV, reaction_prob_VS_etnod_precision_plot, ncol=2, rel_widths = c(4,1))

date = format(Sys.time(), "%a_%b_%d_%X_%Y")
cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=600, width = 250, height = 250, units = 'mm', scale=1)
cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.pdf', sep='', collapse=''), publication_plot, limitsize = F, dpi=600, width = 250, height = 250, units = 'mm', scale=1)

date = format(Sys.time(), "%a_%b_%d_%X_%Y")
cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=200, width = 250, height = 250, units = 'mm', scale=1)

date = format(Sys.time(), "%a_%b_%d_%X_%Y")
cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=100, width = 250, height = 250, units = 'mm', scale=1)
