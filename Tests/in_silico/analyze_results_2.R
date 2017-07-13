library(tidyverse)

D = read.csv('wynyky.matteo') %>% tbl_df

# TO DO : calculate a simple euclidean distance between the real and estimated probabilities instead of all this shit below.

DD =  D %>% 
      select( mols_no, sigma, ETnoD, PTR, contains('ETnoD'), contains('PTR')) %>%
      select(-contains('count') ) %>%
      select(-contains('frag') ) %>%
      select(-contains('precursor') ) 
            

norm_dist = function(etnod,ptr,etnod_estim,ptr_estim) sqrt( ((etnod-etnod_estim)^2 + (ptr-ptr_estim)^2)/2 )
norm_signed_dist = function(etnod,ptr,etnod_estim,ptr_estim) ifelse(ptr < ptr_estim, 1, -1) * sqrt( ((etnod-etnod_estim)^2 + (ptr-ptr_estim)^2)/2 )

  
DDD = 
  DD %>%
  filter(ETnoD+PTR>0) %>%
  mutate(
    ETnoD = ETnoD/(ETnoD+PTR),
    PTR   = PTR/(ETnoD+PTR),
    Basic        = norm_dist(ETnoD, PTR, basic.prob.ETnoD, basic.prob.PTR),
    Intermediate = norm_dist(ETnoD, PTR, intermediate.prob.ETnoD, intermediate.prob.PTR),
    Advanced     = norm_dist(ETnoD, PTR, advanced.prob.ETnoD, advanced.prob.PTR),
    Basic_signed        = norm_signed_dist(ETnoD, PTR, basic.prob.ETnoD, basic.prob.PTR),
    Intermediate_signed = norm_signed_dist(ETnoD, PTR, intermediate.prob.ETnoD, intermediate.prob.PTR),
    Advanced_signed     = norm_signed_dist(ETnoD, PTR, advanced.prob.ETnoD, advanced.prob.PTR)
  ) 

DDD_unsigned = 
  DDD %>%
  select( mols_no, sigma, Basic, Intermediate, Advanced) %>%
  rename( Molecules = mols_no ) %>%
  gather( 'Algorithm', 'dist', -c(1,2) ) 
  

prob_error_plot =
  DDD_unsigned %>%
  ggplot( aes(x=factor(round(sigma,4)), y=dist, color=Algorithm ) ) +
  geom_boxplot() +
  scale_y_continuous(
        limits = c(-.02,1.02),
        labels = scales::percent,
        breaks = (0:10)/10
      ) +
  facet_grid( .~Molecules, labeller = label_both ) +
  theme_minimal() +
  xlab('Standard Deviation of Mass Imprecision') +
  ylab('Estimation Error')
  
cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/prob_error_plot.pdf', prob_error_plot, limitsize = F, dpi=600, width = 250, height = 125, units = 'mm', scale=1)
cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/prob_error_plot.png', prob_error_plot, limitsize = F, dpi=600, width = 250, height = 125, units = 'mm', scale=1)



DDD_signed = 
  DDD %>%
  select( mols_no, sigma, Basic_signed, Intermediate_signed, Advanced_signed) %>%
  rename( Basic = Basic_signed, Intermediate = Intermediate_signed, Advanced = Advanced_signed, Molecules = mols_no ) %>%
  gather( 'Algorithm', 'dist', -c(1,2) ) 

signed_prob_error_plot = 
  DDD_signed %>%
  ggplot( aes(x=factor(round(sigma,4)), y=dist, color=Algorithm ) ) +
  geom_hline( yintercept = 0, size = 1, color='grey' ) +
  geom_boxplot() +
  scale_y_continuous(
    limits = c(-1.02,1.02),
    labels = scales::percent,
    breaks = (-5:5)/5
  ) +
  facet_grid( .~Molecules, labeller = label_both ) +
  theme_minimal() +
  xlab('Standard Deviation of Mass Imprecision') +
  ylab('Estimation Error')

cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/signed_prob_error_plot.pdf', signed_prob_error_plot, limitsize = F, dpi=600, width = 250, height = 125, units = 'mm', scale=1)
cowplot::ggsave('/Users/matteo/Documents/MassTodon/Paper/images/signed_prob_error_plot.png', signed_prob_error_plot, limitsize = F, dpi=600, width = 250, height = 125, units = 'mm', scale=1)
