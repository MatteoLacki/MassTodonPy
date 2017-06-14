library(tidyverse)
library(ggthemes)
library(stringr)

boostrap_res_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_CSV/'
real_res_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_REAL/real.csv'

D = list.files(boostrap_res_path) %>% lapply(read_csv) %>% lapply(tbl_df) %>% bind_rows() %>% mutate(status='bootstrap')
R = read_csv(real_res_path) %>% tbl_df() %>% mutate(status = 'real')

ALL = bind_rows(D,R)
find_among_colnames = function(what, data) colnames(data) %>% .[str_detect(.,what)]

prob_columns  = find_among_colnames('prob', ALL)
count_columns = find_among_colnames('count', ALL)
anion_no_cation_columns = find_among_colnames('anion_did_not_approach_cation', ALL)
info_columns  = c('status', 'ID','WH','WV')
error_columns = c('overestimates','underestimates',
                 'L1_error/original_total_intensity',
                 'L2_error/original_total_intensity',
                 "underestimates/trimmed_total_intensity",
                 "overestimates/trimmed_total_intensity")

DD = ALL[,c(info_columns, error_columns, prob_columns)] 

WV_300 = DD %>% filter(WV==300)

ETnoD_cols = find_among_colnames('ETnoD', WV_300)

WV_300_ETnoD = WV_300[,c(info_columns, ETnoD_cols)]
WV_300_Errors= WV_300[,c(info_columns, error_columns)]
WV_300_anion_no_cation = WV_300[,c(info_columns, anion_no_cation_columns)]

get_algo = Vectorize(function(string) strsplit(string,'_')[[1]][1]) 
check_if_precursor = Vectorize(function(string) tail(strsplit(string,'_')[[1]],1) == 'precursor') 

WV_300_ETnoD = 
  WV_300_ETnoD %>% 
  gather_('algo', 'val', ETnoD_cols) %>%
  filter(!check_if_precursor(algo))  %>%
  mutate(algo = get_algo(algo))

WV_300_ETnoD = WV_300_ETnoD %>% split(.$status)

WV_300_ETnoD$real = 
  WV_300_ETnoD$real %>% group_by(ID, WH, WV) %>% summarize(val=first(val))

WV_300_ETnoD_plot = 
  WV_300_ETnoD$bootstrap %>%
  ggplot(aes(x=factor(WH), y=val, color=algo)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = WV_300_ETnoD$real, color = 'red', shape=19) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Wave Height (Wave Velocity set to 300)') + 
  ylab('Probability of ETnoD')


WV_300_anion_no_cation =  
  WV_300_anion_no_cation %>%
  gather_('algo', 'val', anion_no_cation_columns) %>%
  mutate(algo = get_algo(algo)) %>% split(.$status)


# WV_300_anion_no_cation$bootstrap %>%
#   ggplot(aes(x=factor(WH), y=1-val)) +
#   theme_minimal() +
#   geom_boxplot() +
#   geom_point(data = WV_300_anion_no_cation$real, shape=19) +
#   scale_y_continuous(labels=scales::percent)+
#   xlab('Wave Height (Wave Velocity set to 300)') +
#   ylab('Probability of Fragmentation') +
#   facet_grid(algo~.)

WV_300_anion_no_cation_plot = 
  WV_300_anion_no_cation$real %>%
  ggplot(aes(x=factor(WH), y=1-val, color=algo)) +
  theme_minimal() +
  geom_point(shape=19) +
  scale_y_continuous(labels=scales::percent)+
  xlab('Wave Height (Wave Velocity set to 300)') +
  ylab('Anion meets Cation') + 
  theme(legend.position="bottom")

cowplot::plot_grid(WV_300_ETnoD_plot, WV_300_anion_no_cation_plot, nrow=2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

WH_150 = DD %>% filter(WH==150)

ETnoD_cols = find_among_colnames('ETnoD', WH_150)

WH_150_ETnoD = WH_150[,c(info_columns, ETnoD_cols)]
WH_150_Errors= WH_150[,c(info_columns, error_columns)]
WH_150_NoFrags = WH_150[,c(info_columns, nofrag_columns)]

get_algo = Vectorize(function(string) strsplit(string,'_')[[1]][1]) 
check_if_precursor = Vectorize(function(string) tail(strsplit(string,'_')[[1]],1) == 'precursor') 

WH_150_ETnoD = 
  WH_150_ETnoD %>% 
  gather_('algo', 'val', ETnoD_cols) %>%
  filter(!check_if_precursor(algo))  %>%
  mutate(algo = get_algo(algo))

WH_150_ETnoD = WH_150_ETnoD %>% split(.$status)

WH_150_ETnoD$bootstrap %>%
  ggplot(aes(x=factor(WV), y=val)) + 
  theme_minimal() +
  geom_boxplot() +
  geom_point(data = WH_150_ETnoD$real, color = 'red', shape=19) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(labels=scales::percent)+
  xlab('Wave Velocity (Wave Height set to 150)') + 
  ylab('Probability of ETnoD')




# ggplot(aes(x=X, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
#   geom_boxplot() +
#   theme_minimal() +
#   xlab(xlab_label)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position = 'top' ) +
#   scale_y_continuous(labels=scales::percent) 
# 
# 
# 
# make_box_plot = function(D, xlab_label, x_choice, x) D %>% 
#   filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
#   select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
#   mutate(WH = factor(WH), WV = factor(WV) ) %>%
#   mutate_(X=x) %>%
#   ggplot(aes(x=X, y=ETnoD, fill=algo, alpha =anion_approached_cation)) +
#   geom_boxplot() +
#   theme_minimal() +
#   xlab(xlab_label)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         legend.position = 'top' ) +
#   scale_y_continuous(labels=scales::percent) 
# 
# make_react_plot = function(D, xlab_label, x_choice, x) D %>% 
#   filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
#   select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
#   mutate(WH = factor(WH), WV = factor(WV) ) %>%
#   mutate_(X=x) %>%
#   ggplot(aes(x=X, y=anion_approached_cation)) +
#   geom_boxplot()+
#   theme_minimal() +
#   xlab(xlab_label)+
#   ylab('Probability of Anion Reaction') +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_y_continuous(labels=scales::percent) 
#   
# quant = function(vals, prob) 
#   if(any(is.na(vals))) return(NA) else return(quantile(vals, probs=prob)) 
# 
# make_error_plot = function(D, xlab_label, x_choice, x, min_prob, max_prob){
#   LD= D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'") %>%
#     select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
#     mutate(WH = factor(WH), WV = factor(WV) ) %>%
#     mutate_(X=x) %>%
#     group_by(X, algo) %>%
#     summarize(mean_val = mean(ETnoD),
#               prob_max = quant(ETnoD, max_prob),
#               prob_min = quant(ETnoD, min_prob) ) %>% ungroup
#   RD = D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='real'") %>%
#     select(WH, WV, algo, ETnoD, PTR, ID, anion_approached_cation) %>%
#     mutate(WH = factor(WH), WV = factor(WV) ) %>%
#     mutate_(X=x) %>%
#     group_by(X, algo) 
#   ggplot(LD) +
#     geom_crossbar(
#       aes(x = X,
#           y = mean_val,
#           ymin = prob_min, 
#           ymax = prob_max,
#           color = algo), 
#       width = 0.4,
#       position = position_dodge()
#     )+
#     geom_point(
#       data = RD,
#       aes(
#         x = X,
#         y = ETnoD
#       )
#     ) +
#     theme_minimal() +
#     xlab(xlab_label)+
#     ylab('ETnoD Bootstrap Estimates')+
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 1),
#       legend.position = 'top'
#     ) +
#     scale_y_continuous(labels=scales::percent)  }
# 
# make_plot = function(D, x_arm="WV", min_prob, max_prob ){
#   if( x_arm=="WV" ){ 
#           xlab_label = 'Wave Height (Wave Velocity set to 300)'
#           x_choice = "WV==300"
#           x = 'WH'
#   } else {xlab_label = 'Wave Velocity (Wave Height set to 150)'
#           x_choice = "WH==150"
#           x = 'WV' }
#   
#   PTR_ETnoD  = make_box_plot(D, xlab_label, x_choice, x)
#   prob_react = make_react_plot(D, xlab_label, x_choice, x)
#   error_plot = make_error_plot(D, xlab_label, x_choice, x, min_prob, max_prob)
#   
#   list( 
#     # PTR_ETnoD = PTR_ETnoD, 
#     error     = error_plot, 
#     prob_react= prob_react 
#   )
# }
# 
# min_prob = .025 
# max_prob = .975
# 
# plot_grid(plotlist=make_plot(D, "WV", min_prob, max_prob), nrow=2)
# plot_grid(plotlist=make_plot(D, "WH", min_prob, max_prob), nrow=2)
# 
# 
# 
# make_plot_frag = function(D, x_arm="WV", fasta='RPKPQQFFGLM' ){
#   if( x_arm=="WV" ){ 
#     xlab_label = 'Wave Height (Wave Velocity set to 300)'
#     x_choice = "WV==300"
#     x = 'WH'
#   } else {xlab_label = 'Wave Velocity (Wave Height set to 150)'
#   x_choice = "WH==150"
#   x = 'WV' }
# 
#   PD = D %>% filter_(x_choice, "count_or_prob=='prob'", "real_or_sim=='sim'", "algo=='base'") %>%
#         select( -PTR, -PTR_precursor, -anion_approached_cation, -ETnoD_precursor, -ETnoD, -algo,
#                 -anion_did_not_approach_cation, -fragmentation, -no.fragmentation, -reactions, 
#                 -real_or_sim, -total_frags, - total_reactions, -unreacted_precursors) %>%
#         gather('frag_place', 'frag_prob', 1:9) %>%
#         filter_(x_choice) 
# 
#   PD_dummy = data.frame(ID=NA, WH=NA, WV=NA, count_or_prob=NA, frag_place=paste0('X',0:10, sep=''), frag_prob=NA )  
#   PD = bind_rows(PD, PD_dummy)
#   PD = PD %>% mutate( WH = factor(WH), 
#                       WV = factor(WV),
#                       frag_place = ordered(frag_place, levels=paste0('X',0:10, sep='') ) 
#               ) %>%
#               mutate_(X=x) 
#   
#   
#   plot_AA <- Vectorize(function(x) strsplit(fasta,'')[[1]][as.integer(x)+1])
#  
#   look_up_table = plot_AA(0:10)
#   names(look_up_table) = paste0('X',0:10)
#   
#   PD %>%
#     filter(!is.na(X)) %>%
#     ggplot(aes(x=X,y=frag_prob)) +
#     geom_boxplot() + 
#     facet_grid( .~frag_place, drop = FALSE, labeller=as_labeller(look_up_table) ) +
#     theme_minimal()+
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),
#           legend.position = 'top' ) + 
#     scale_y_continuous(labels=scales::percent) + 
#     ylab('Probability of Fragmentation')+
#     xlab(xlab_label)
# }
# 
# make_plot_frag(D, x_arm="WV", fasta='RPKPQQFFGLM' )
# make_plot_frag(D, x_arm="WH", fasta='RPKPQQFFGLM' )
# 
# 
