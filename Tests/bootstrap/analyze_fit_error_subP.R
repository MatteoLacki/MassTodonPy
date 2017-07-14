library('tidyverse'); library('stringr'); library('car')

save_plot = function(a_particular_plot, name, dpi=600, width=300, height=300){
  cowplot::ggsave(
    paste0('/Users/matteo/Documents/MassTodon/Paper/images/',name, sep='', collapse=''),
    a_particular_plot,
    limitsize = F,
    dpi=dpi, width =width, height = height, units = 'mm', scale=1)
  print('Saved')
}

D = list.files('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/CSV_03_07_2017_mzPrec-065', full.names = T) %>%
  lapply(read.csv) %>%
  lapply(as_data_frame) %>%
  bind_rows %>%
  mutate(WH = factor(WH), WV = factor(WV) )

D150 = D %>% filter(WH == '150')
D300 = D %>% filter(WV == '300')

D150 = D150 %>% select(ID, WV, real_or_bootstrap, L1_error.intensity_within_tolerance, L1_fit_error_and_unused_intensity.original_total_intensity)
D300 = D300 %>% select(ID, WH, real_or_bootstrap, L1_error.intensity_within_tolerance, L1_fit_error_and_unused_intensity.original_total_intensity)

# Binary plot search
plot_150_v =
  D150 %>%
  gather('error_type', 'L1_error', 4:5) %>%
  mutate(
    error_type = plyr::revalue( error_type,
                           c( 'L1_error.intensity_within_tolerance'='fit error within scooped area',
                              'L1_fit_error_and_unused_intensity.original_total_intensity' = 'mismatch error'))
  ) %>%
  ggplot( aes( x=WV, y=L1_error, color=error_type ) ) +
  geom_boxplot() +
  geom_violin() +
  coord_flip() +
  theme_minimal()+
  theme(legend.position ="top")+
  scale_y_continuous(
    labels = scales::percent, position = "right",
    limits = c(0,1),,
    breaks = (0:10)/10
  ) +
  xlab('Wave Velocity (WH = 1,5)') +
  ylab('Relative Error')+
  geom_hline(yintercept = 0, color='grey', size=2, alpha=.5) +
  geom_hline(yintercept = 1, color='grey', size=2, alpha=.5) +
  guides(color=guide_legend(title=""))


plot_300_v =
  D300 %>%
  gather('error_type', 'L1_error', 4:5) %>%
  mutate(
    error_type = plyr::revalue( error_type,
                                c( 'L1_error.intensity_within_tolerance'='fit error within scooped area',
                                   'L1_fit_error_and_unused_intensity.original_total_intensity' = 'mismatch error'))
  ) %>%
  ggplot( aes( x=WH, y=L1_error, color=error_type ) ) +
  geom_boxplot() +
  geom_violin() +
  coord_flip() +
  theme_minimal()+
  theme(legend.position ="top")+
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0,1),
    breaks = (0:10)/10
  ) +
  geom_hline(yintercept = 0, color='grey', size=2, alpha=.5) +
  geom_hline(yintercept = 1, color='grey', size=2, alpha=.5) +
  xlab('Wave Height (WV = 300)') +
  ylab('Relative Error')+
  guides(color=FALSE)

big_plot = cowplot::plot_grid(plot_150_v, plot_300_v, nrow=2, align='v', rel_heights = c(2,.8) )
save_plot(big_plot, 'fit_error.png', width=150)
save_plot(big_plot, 'fit_error.pdf', width=150)
