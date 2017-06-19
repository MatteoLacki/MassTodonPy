library(tidyverse)
library(ggthemes)
library(stringr)

boostrap_res_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_CSV/'
real_res_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_REAL/real.csv'


D = list.files(boostrap_res_path, full.names = T ) %>% lapply(read_csv) %>% lapply(tbl_df) %>% bind_rows() %>% mutate(status='bootstrap')
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

total_reaction_columns = 'inter_count_total_reactions'

DD = ALL[,c(info_columns, error_columns, prob_columns)]

