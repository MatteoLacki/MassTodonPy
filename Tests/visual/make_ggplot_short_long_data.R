library(tidyverse); library(ggthemes); library(ggiraph); library(ggthemes); library(htmlwidgets)

D = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_long.csv')
S = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_short.csv') 

all((    S %>% 
             mutate( mz = (mz_L+mz_R)/2) %>%
             group_by(mz) %>%
             summarize( 
                 tot_intensity = sum(tot_intensity),
                 n = n()
             ) %>% 
             filter( tot_intensity > 0 ) %>% 
             data.frame() %>%
             select(n)
) == 1)

S = S %>% 
    mutate( mz = (mz_L+mz_R)/2) %>%
    group_by(mz) %>%
    summarize( 
        mz_L = first(mz_L),
        mz_R = first(mz_R),
        tot_intensity = sum(tot_intensity),
        tot_estimate  = sum(tot_estimate)
    ) # Things with the same mz are aggregated: these are only artificially added 0 intensity G nodes.



min_mz = as.integer(min(D$mz_L, na.rm = T)) - 1
max_mz = as.integer(max(D$mz_R, na.rm = T)) + 1

# base_plot = 
#     D %>%    
#     drop_na() %>%
#     ggplot() +
#     geom_rect( aes(xmin = mz_L, xmax = (mz_R+mz_L)/2, ymin = 0, ymax = tot_intensity ), fill = 'white' ) +
#     geom_rect( aes(xmin = (mz_R+mz_L)/2, xmax = mz_R, ymin = 0, ymax = tot_estimate ),  fill = 'red' ) +
#     theme_hc(bgcolor = "darkunica") +
#     scale_x_continuous(
#         breaks = seq( min_mz, max_mz, by = 2 )
#     ) +
#     xlab('m/z')+
#     ylab('intensity')
# 
# ggsave( filename = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_plot.pdf',
#         plot     = base_plot,
#         limitsize= F,
#         height   = 10,
#         width    = 1500  )
