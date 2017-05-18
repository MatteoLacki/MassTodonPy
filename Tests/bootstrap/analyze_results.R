library(tidyverse)

R = read.csv('real_data.csv') %>% tbl_df()
S = read.csv('simulations.csv') %>% tbl_df()

D = bind_rows(R,S) 
WV_300_D = D %>% select(WH, WV, algo, ETnoD, PTR, real_or_sim, ID) %>% filter(WV==300)
WV_300_D %>% filter(real_or_sim=='real', WH==150, WV==300)



WV_300_plot = S %>% select(WH, WV, algo, ETnoD, PTR, ID) %>% filter(WV==300) %>%
  mutate(WH = factor(WH)) %>%
  ggplot(aes(x=WH, y=ETnoD, color=algo)) +
  geom_boxplot() +
  theme_minimal() +
  xlab('Wave Height (Wave Velocity set to 300)')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()

  
  
# WH_150_plot = 
  S %>% select(WH, WV, algo, ETnoD, PTR, ID) %>% filter(WH==150) %>%
  mutate(WV = factor(WV)) %>%
  ggplot(aes(x=WV, y=ETnoD, color=algo)) +
  geom_boxplot() +
  theme_minimal() +
  xlab('Wave Velocity (Wave Height set to 150)')+
  scale_y_continuous(labels=scales::percent) +
  coord_flip()

library(cowplot)

plot_grid( WV_300_plot, WH_150_plot)
