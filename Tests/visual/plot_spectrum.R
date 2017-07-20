library(tidyverse)
library(ggthemes)
library(jsonlite)

all_ubis = 
    list.files('ubiquitins', full.names = T) %>%
    lapply( read_csv )




OD = read_json('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_original.json') 
OD = data_frame(mz        = OD %>% sapply('[[',1),
                intensity = OD %>% sapply('[[',2) )



OD %>% plot(type='h')

ubi_data = read_csv('ubi_data.csv')




# 
# D = ubi_data %>% split(.$where)
# 
# D$not_explainable %>% 
#     mutate(mz = (mz_L + mz_R)/2) %>%
#     select(mz, tot_intensity ) %>%
#     plot(tyoe='h')
# 
# left = 1000; right = 2160
# u_lr = ubi_data %>% filter(mz_L > left, mz_L < right)
# 
# 
# ubi_data %>% 
#     mutate( mz = (mz_L + mz_R)/2) %>%
#     select( mz, tot_intensity ) %>% plot(type='l')
# 
# ubi_data %>%
#     mutate( tot_estimate = -tot_estimate ) %>%
#     gather( 'tot', 'value', 9:10 ) %>%
#     ggplot( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = value, fill = tot)) +
#     geom_rect() + 
#     geom_rangeframe() +
#     theme_tufte() 
# 
# 
# 
# 
# draw_fragment = function(left, right, data){
#     D = data %>% filter(mz_L > left, mz_L < right) 
#     
#     plot 
#         ggplot( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = estimate), fill = 'red', color = 'red' ) +
#         geom_rect() + 
#         geom_rangeframe() +
#         theme_tufte() 
#     
# }
# 
# draw_fragment(2135, 2160, D$theoretically_achievable)
# 
# 
# D$theoretically_achievable %>% 
#     ggplot( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = estimate), fill = 'white' ) +
#     geom_rect() + 
#     theme_dark()
# 
