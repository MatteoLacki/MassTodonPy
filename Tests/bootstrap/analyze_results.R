library(tidyverse)
library(cowplot)

S = read.csv('simulations_3.csv') %>% tbl_df()
R = read.csv('real_data_3.csv') %>% tbl_df()

D = bind_rows(R,S) 

WV_300_D = D$prob %>% select(WH, WV, algo, ETnoD, PTR, real_or_sim, ID) %>% filter(WV==300)
WV_300_D %>% filter(real_or_sim=='real', WH==150, WV==300)

make_plot = function(D, x_arm="WV" ){
  if( x_arm=="WV" ){ 
          xlab_label = 'Wave Height (Wave Velocity set to 300)'
          x_choice = "WV==300"
          x = 'WH'
  } else {xlab_label = 'Wave Velocity (Wave Height set to 150)'
          x_choice = "WH==150"
          x = 'WV' }
  PTR_ETnoD = D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
              select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
              mutate(WH = factor(WH), WV = factor(WV) ) %>%
              mutate_(X=x) %>%
              ggplot(aes(x=X, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
              geom_boxplot() +
              theme_minimal() +
              xlab(xlab_label)+
              theme(
                axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = 'top'
              ) +
              scale_y_continuous(labels=scales::percent) 
  
  prob_react= D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
              select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
              mutate(WH = factor(WH), WV = factor(WV) ) %>%
              mutate_(X=x) %>%
              ggplot(aes(x=X, y=anion_approached_cation)) +
              geom_boxplot()+
              theme_minimal() +
              xlab(xlab_label)+
              theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
              scale_y_continuous(labels=scales::percent) 
  
  list(PTR_ETnoD=PTR_ETnoD, prob_react=prob_react)
}

plot_grid(plotlist=make_plot(D, "WV"), nrow=2)
plot_grid(plotlist=make_plot(D, "WH"), nrow=2)


x_arm="WV"
