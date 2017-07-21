library(tidyverse)
library(ggthemes)
library(jsonlite)

all_ubis = 
    list.files('ubiquitins', full.names = T) %>%
    lapply( read_csv )

spectrum = all_ubis[[100]]

filenames = list.files('ubiquitins') %>% 
    sapply(tools::file_path_sans_ext) %>% 
    sapply(function(x) paste0( '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiquitins_plots/', x, '.pdf', sep='', collapse='') )

for( i in 1:length(all_ubis)){
    spectrum = all_ubis[[i]]
    pdf(filenames[i])
        plot(spectrum$mass, spectrum$intensity, type='h' )
    dev.off()
}

# plot_stupid = 
#     spectrum %>%
#     ggplot( aes(x=mass, xend=mass, y=0, yend=intensity )) +
#     geom_segment() + 
#     theme_minimal() + 
#     theme(
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#     ) +
#     xlab('m/z') +
#     ylab('intensity') + 
#     scale_x_continuous(
#         breaks = seq( min_mz, max_mz, by = 1 )
#     ) 
#     
# ggsave( filename = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiquitins_plots/test.pdf', 
#         limitsize = F,
#         width = 5000,
#         height = 30,
#         plot_stupid)
# 
# 
# 
OD = read_json('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_original.json')
OD = data_frame(mz        = OD %>% sapply('[[',1),
                intensity = OD %>% sapply('[[',2) )

plot(OD, type='h')

hchart(OD, "scatter", x = day_of_week, y = diff_ppt)
# 
# 
# 
# OD %>% plot(type='h')
# 
# ubi_data = read_csv('ubi_data.csv')




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
