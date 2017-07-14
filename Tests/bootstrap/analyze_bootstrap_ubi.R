library(tidyverse); library(stringr); library(car); library(jsonlite)


experimentalParameters %>% data.frame
initial_run_files = list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun') %>% tools::file_path_sans_ext()
ID_files = 
  read_json('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ID_file.json') %>% 
  lapply(function(x) data.frame(file=tools::file_path_sans_ext(x)) ) %>%
  bind_rows(.id = 'ID') %>% 
  mutate(ID = as.integer(ID)) %>%
  tbl_df()
initial_run_ID = (ID_files %>% filter( file %in% initial_run_files ))$ID

D = list.files('UBI_csv', full.names = T) %>%
  lapply(read.csv) %>%
  lapply(as_data_frame) 

E = lapply(D, function(d) d %>% mutate(preActive = factor(preActive))) %>% bind_rows(.id='ID_2') 


E %>% filter( preActive == 'none' ) %>%
  group_by( supActive, preActive, retentionTime,precursorMZ  ) %>%
  summarize( n = n() ) %>% data.frame
## check if IDs within the boundries

E %>% filter( 
  preActive == 'none', 
  precursorMZ %in% c(952, 1428),
  ! ID %in% initial_run_ID
) %>%
group_by( supActive, preActive, retentionTime,precursorMZ  ) %>%
summarize( n = n() ) %>% 
filter(n>251)  

#####
E %>% filter( preActive != 'none' ) %>%
  group_by( supActive, preActive, retentionTime ) %>%
  summarize( n = n() ) %>% data.frame %>% 
  filter(n>251)  


(E %>% filter( 
  preActive == 'none', 
  precursorMZ %in% c(952, 1428),
  ! ID %in% initial_run_ID
) %>% nrow) +
(E %>% filter( preActive != 'none' ) %>% nrow)



################## FINAL DATA PREPARATION
DE = E %>% 
  filter(
    precursorMZ %in% c(952, 1428),
    ! ID %in% initial_run_ID
  )


check_nrow = function(DE) DE %>% group_by( supActive, preActive, retentionTime, precursorMZ ) %>% summarize( n = n() ) %>% data.frame %>% filter(n>251) %>% nrow
assertthat::are_equal(DE %>% group_by( supActive, preActive, retentionTime, precursorMZ ) %>% summarize( n = n() ) %>% data.frame %>% filter(n>251) %>% nrow, 0)


# D3 = D2 %>% mutate(
#   ID = factor(ID), 
#   preActive = ordered(preActive, levels = c('none', '15', '20')),
#   supActive = factor(supActive, levels = c('Supplementary Activation OFF', 'Supplementary Activation ON'), labels = c('ON', 'OFF') ),
#   Q = as.integer(8500/precursorMZ)
# )

# D3 = D2 %>%
#   filter(precursorMZ %in% c(952, 1428)) %>%
#   mutate(Q = as.integer(8500/precursorMZ)) %>%
#   split(.$Q)

DE = DE %>% mutate(Q = as.integer(8500/precursorMZ)) %>% split(.$Q)


DD8 =
  DE$'8' %>%
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
  split(.$real_or_bootstrap)

estimates_plot_8 =
  DD8$boot %>%
  ggplot(aes(x=ordered(retentionTime), y=prob_val, fill = Algorithm)) +
  theme_grey() +
  geom_boxplot() +
  geom_point(data = DD8$real, color = 'red', shape=19) +
  facet_grid(supActive ~ preActive, labeller = label_both) +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-.,
      labels= scales::percent,
      name  = "PTR"
    )
  ) +
  theme(legend.position="top")+
  xlab('Retention Time') +
  ylab('ETnoD')



DD5 =
  DE$'5' %>%
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
  split(.$real_or_bootstrap)

estimates_plot_5 =
  DD5$boot %>%
  ggplot(aes(x=ordered(retentionTime), y=prob_val, fill = Algorithm)) +
  theme_grey() +
  geom_boxplot() +
  geom_point(data = DD5$real, color = 'red', shape=19) +
  facet_grid(supActive ~ preActive, labeller = label_both) +
  scale_y_continuous(
    labels=scales::percent,
    sec.axis = sec_axis(
      trans = ~1-.,
      labels= scales::percent,
      name  = "PTR"
    )
  ) +
  theme(legend.position="bottom")+
  xlab('Retention Time') +
  ylab('ETnoD')

estimates_plot = cowplot::plot_grid(estimates_plot_8, estimates_plot_5, nrow=2, align='v', labels = c('Q = 8', 'Q = 5'))

# DD8$boot %>%
#   group_by( supActive, preActive, retentionTime, Algorithm ) %>%
#   summarize( n = n() ) %>% data.frame


# 
# 
# 
# 
# 
# capFirst <- function(s) {
#   paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
# }
# get_algo = Vectorize(function(x) capFirst(strsplit(x,'[.]')[[1]][1]))
# 
# reaction_D150 = D150$boot %>%
#   gather('algo','prob_reaction',1:3) %>%
#   mutate(Algorithm = get_algo(algo))
# 
# reaction_plot150 =
#   reaction_D150 %>%
#   ggplot(aes(x=WV, y=prob_reaction, fill=Algorithm, color=Algorithm)) +
#   geom_boxplot()+
#   theme_minimal() +
#   scale_y_continuous(
#     labels=scales::percent,
#     sec.axis = sec_axis(
#       trans = ~1-.,
#       labels=scales::percent,
#       name  = "No Reaction"
#     )
#   )+
#   scale_x_discrete(position = "top") +
#   xlab('Wave Velocity (Wave Height set to 1,5)') +
#   ylab('Reaction') +
#   theme(legend.position="top")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# WH150 = cowplot::plot_grid(reaction_plot150, estimates_plot150, nrow=2, align='v')
# 
# ################################################################################################
# 
# D300 = D %>% filter(WV == '300') %>%
#   select( contains('prob.ETnoD'), contains('anion_approached_cation'), real_or_bootstrap, WH, WV, ID, intensity_within_tolerance ) %>%
#   gather("prob_type", "prob_val", 1:6 ) %>%
#   split(.$real_or_bootstrap)
# 
# 
# estimates_plot300 =
#   D300$boot %>%
#   ggplot(aes(x=WH, y=prob_val)) +
#   theme_minimal() +
#   geom_boxplot() +
#   geom_point(data = D300$real, color = 'red', shape=19) +
#   theme( axis.text.x = element_text(angle = 90, hjust = 1)
#          # axis.title.x=element_blank()
#          # axis.text.x=element_blank()
#   ) +
#   scale_y_continuous(
#     labels   = scales::percent,
#     sec.axis = sec_axis(
#       trans = ~1-.,
#       labels=scales::percent,
#       name  = 'PTR'
#     )
#   )+
#   xlab('Wave Height (Wave Velocity set to 300)') +
#   ylab('ETnoD')
# 
# reaction_D300 = D300$boot %>%
#   gather('algo','prob_reaction',1:3) %>%
#   mutate(Algorithm = get_algo(algo))
# 
# reaction_plot300 =
#   reaction_D300 %>%
#   ggplot(aes(x=WH, y=prob_reaction, fill=Algorithm, color=Algorithm)) +
#   geom_boxplot()+
#   theme_minimal() +
#   scale_y_continuous(
#     labels=scales::percent,
#     sec.axis = sec_axis(
#       trans = ~1-.,
#       labels=scales::percent,
#       name  = "No Reaction"
#     )
#   )+
#   scale_x_discrete(position = "top") +
#   xlab('Wave Height (Wave Velocity set to 300)') +
#   ylab('Reaction') +
#   theme(legend.position="top")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# WV300 = cowplot::plot_grid(reaction_plot300, estimates_plot300, nrow=2, align='v')
# 
# 
# ################################################################################################
# 
# reaction_prob_VS_etnod_precision_plot =
#   D %>%
#   select( contains('prob.ETnoD'), contains('anion_approached_cation'), real_or_bootstrap, WH, WV, ID, intensity_within_tolerance ) %>%
#   gather("prob_type", "prob_val", 1:6 ) %>%
#   gather("algo", "prob_reaction", 1:3)  %>%
#   drop_na %>% group_by(ID, prob_type, WH) %>%
#   mutate( sd_estim = sd(prob_val),
#           var_estim = var(prob_val),
#           ave_intensity_within_tol = mean(intensity_within_tolerance),
#           ave_anion_approached_cation = mean(prob_reaction)) %>%
#   ungroup %>% select(-real_or_bootstrap, -prob_type, -WH, -WV) %>%
#   mutate(Algorithm = get_algo(algo)) %>%
#   ggplot(aes(x=prob_reaction, y=1.0/var_estim, color=Algorithm, fill=Algorithm)) +
#   geom_point() +
#   # geom_hex(bins=100) +
#   geom_smooth(method='lm', formula=y ~ poly(x, 2)+poly(1/(1-x),4)) +
#   theme_minimal() +
#   scale_x_continuous(
#     labels=scales::percent,
#     breaks=seq(0,1,.05)
#   )+
#   xlab('Average Probability of a Reaction') +
#   ylab('Precision of ETnoD Estimates')+
#   theme(legend.position="bottom")
# 
# 
# ################################################################################################
# 
# WH_WV = cowplot::plot_grid(WH150, WV300, ncol=2, align='h', rel_widths = c(1.5,.75))
# 
# publication_plot = cowplot::plot_grid(WH_WV, reaction_prob_VS_etnod_precision_plot, nrow=2, rel_heights = c(1,.5))
# # publication_plot = cowplot::plot_grid(WH_WV, reaction_prob_VS_etnod_precision_plot, ncol=2, rel_widths = c(4,1))
# 
# date = format(Sys.time(), "%a_%b_%d_%X_%Y")
# cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=600, width = 250, height = 250, units = 'mm', scale=1)
# cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.pdf', sep='', collapse=''), publication_plot, limitsize = F, dpi=600, width = 250, height = 250, units = 'mm', scale=1)
# 
# date = format(Sys.time(), "%a_%b_%d_%X_%Y")
# cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=200, width = 250, height = 250, units = 'mm', scale=1)
# 
# date = format(Sys.time(), "%a_%b_%d_%X_%Y")
# cowplot::ggsave(paste0('ETnoD_PTR_meetingCation___', date, '.png', sep='', collapse=''), publication_plot, limitsize = F, dpi=100, width = 250, height = 250, units = 'mm', scale=1)
