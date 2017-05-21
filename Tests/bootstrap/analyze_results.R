library(tidyverse)

S = read.csv('simulations_3.csv') %>% tbl_df()
R = read.csv('real_data_3.csv') %>% tbl_df()

D = bind_rows(R,S) 
WV_300_D = D %>% select(WH, WV, algo, ETnoD, PTR, real_or_sim, ID) %>% filter(WV==300)
WV_300_D %>% filter(real_or_sim=='real', WH==150, WV==300)


S %>% select(WH, WV, algo, ETnoD, PTR, ID) %>% filter(WV==300)%>% data.frame

WV_300_plot = 
  S %>% 
  select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>% filter(WV==300) %>%
  mutate(WH = factor(WH)) %>%
  ggplot(aes(x=WH, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
  geom_boxplot() +
  theme_minimal() +
  xlab('Wave Height (Wave Velocity set to 300)')+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = 'top'
  ) +
  scale_y_continuous(labels=scales::percent) 
  # coord_flip()

WV_300_plot_anions = 
  S %>% select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>% filter(WV==300) %>%
  mutate(WH = factor(WH)) %>%
  ggplot(aes(x=WH, y=anion_approached_cation)) +
  geom_boxplot()+
  theme_minimal() +
  xlab('Wave Height (Wave Velocity set to 300)')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(labels=scales::percent) 

library(cowplot)
  
plot_grid(WV_300_plot, WV_300_plot_anions, nrow=2)
  
  
WH_150_plot =
  S %>% select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>% filter(WH==150) %>%
  mutate(WV = factor(WV)) %>%
  ggplot(aes(x=WV, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
  geom_boxplot() +
  theme_minimal() +
  xlab('Wave Velocity (Wave Height set to 150)')+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = 'top'
  ) +
  scale_y_continuous(labels=scales::percent) 

WH_150_plot_anions = 
  S %>% select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>% filter(WH==150) %>%
  mutate(WV = factor(WV)) %>%
  ggplot(aes(x=WV, y=anion_approached_cation)) +
  geom_boxplot()+
  theme_minimal() +
  xlab('Wave Velocity (Wave Height set to 150)')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(labels=scales::percent) 


library(cowplot)

plot_grid( WH_150_plot, WH_150_plot_anions, nrow=2)
