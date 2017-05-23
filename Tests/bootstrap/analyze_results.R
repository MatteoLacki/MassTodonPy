library(tidyverse)
library(cowplot)
library(ggthemes)

S = read.csv('simulations_3.csv') %>% tbl_df()
R = read.csv('real_data_3.csv') %>% tbl_df()

D = bind_rows(R,S) 

WV_300_D = D$prob %>% select(WH, WV, algo, ETnoD, PTR, real_or_sim, ID) %>% filter(WV==300)
WV_300_D %>% filter(real_or_sim=='real', WH==150, WV==300)

make_box_plot = function(D, xlab_label, x_choice, x) D %>% 
  filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
  select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
  mutate(WH = factor(WH), WV = factor(WV) ) %>%
  mutate_(X=x) %>%
  ggplot(aes(x=X, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
  geom_boxplot() +
  theme_minimal() +
  xlab(xlab_label)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = 'top' ) +
  scale_y_continuous(labels=scales::percent) 

make_react_plot = function(D, xlab_label, x_choice, x) D %>% 
  filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
  select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
  mutate(WH = factor(WH), WV = factor(WV) ) %>%
  mutate_(X=x) %>%
  ggplot(aes(x=X, y=anion_approached_cation)) +
  geom_boxplot()+
  theme_minimal() +
  xlab(xlab_label)+
  ylab('Probability of Anion Reaction') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(labels=scales::percent) 
  
quant = function(vals, prob) 
  if(any(is.na(vals))) return(NA) else return(quantile(vals, probs=prob)) 

make_error_plot = function(D, xlab_label, x_choice, x, min_prob, max_prob){
  LD= D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
    select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
    mutate(WH = factor(WH), WV = factor(WV) ) %>%
    mutate_(X=x) %>%
    group_by(X, algo) %>%
    summarize(mean_val = mean(ETnoD),
              prob_max = quant(ETnoD, max_prob),
              prob_min = quant(ETnoD, min_prob) ) %>% ungroup
  RD = D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='real'") %>%
    select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
    mutate(WH = factor(WH), WV = factor(WV) ) %>%
    mutate_(X=x) %>%
    group_by(X, algo) 
  ggplot(LD) +
    geom_crossbar(
      aes(x = X,
          y = mean_val,
          ymin = prob_min, 
          ymax = prob_max,
          color = algo), 
      width = 0.4,
      position = position_dodge()
    )+
    geom_point(
      data = RD,
      aes(
        x = X,
        y = ETnoD
      )
    ) +
    theme_minimal() +
    xlab(xlab_label)+
    ylab('ETnoD Bootstrap Estimates')+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = 'top'
    ) +
    scale_y_continuous(labels=scales::percent)  }

make_plot = function(D, x_arm="WV", min_prob, max_prob ){
  if( x_arm=="WV" ){ 
          xlab_label = 'Wave Height (Wave Velocity set to 300)'
          x_choice = "WV==300"
          x = 'WH'
  } else {xlab_label = 'Wave Velocity (Wave Height set to 150)'
          x_choice = "WH==150"
          x = 'WV' }
  
  PTR_ETnoD  = make_box_plot(D, xlab_label, x_choice, x)
  prob_react = make_react_plot(D, xlab_label, x_choice, x)
  error_plot = make_error_plot(D, xlab_label, x_choice, x, min_prob, max_prob)
  
  list( 
    # PTR_ETnoD = PTR_ETnoD, 
    error     = error_plot, 
    prob_react= prob_react 
  )
}

min_prob = .025 
max_prob = .975

plot_grid(plotlist=make_plot(D, "WV", min_prob, max_prob), nrow=2)
plot_grid(plotlist=make_plot(D, "WH", min_prob, max_prob), nrow=2)



make_plot_frag = function(D, x_arm="WV", fasta='RPKPQQFFGLM' ){
  if( x_arm=="WV" ){ 
    xlab_label = 'Wave Height (Wave Velocity set to 300)'
    x_choice = "WV==300"
    x = 'WH'
  } else {xlab_label = 'Wave Velocity (Wave Height set to 150)'
  x_choice = "WH==150"
  x = 'WV' }

  PD = D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'", "algo=='base'") %>%
        select( -PTR, -PTR_precursor, -anion_approached_cation, -ETnoD_precursor, -ETnoD, -algo,
                -anion_did_not_approach_cation, -fragmentation, -no.fragmentation, -reactions, 
                -real_or_sim, -total_frags, - total_reactions, -unreacted_precursors) %>%
        gather('frag_place', 'frag_prob', 1:9) %>%
        filter_(x_choice) 

  PD_dummy = data.frame(ID=NA, WH=NA, WV=NA, count_or_prob=NA, frag_place=paste0('X',0:10, sep=''), frag_prob=NA )  
  PD = bind_rows(PD, PD_dummy)
  PD = PD %>% mutate( WH = factor(WH), 
                      WV = factor(WV),
                      frag_place = ordered(frag_place, levels=paste0('X',0:10, sep='') ) 
              ) %>%
              mutate_(X=x) 
  
  
  plot_AA <- Vectorize(function(x) strsplit(fasta,'')[[1]][as.integer(x)+1])
 
  look_up_table = plot_AA(0:10)
  names(look_up_table) = paste0('X',0:10)
  
  PD %>%
    filter(!is.na(X)) %>%
    ggplot(aes(x=X,y=frag_prob)) +
    geom_boxplot() + 
    facet_grid( .~frag_place, drop = FALSE, labeller=as_labeller(look_up_table) ) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = 'top' ) + 
    scale_y_continuous(labels=scales::percent) + 
    ylab('Probability of Fragmentation')+
    xlab(xlab_label)
}

make_plot_frag(D, x_arm="WV", fasta='RPKPQQFFGLM' )
make_plot_frag(D, x_arm="WH", fasta='RPKPQQFFGLM' )


