library(tidyverse)
library(ggthemes)

D = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data.csv')



min_mz = as.integee(min(D$mz_L)) - 1
min_mz = as.integee(max(D$mz_R)) + 1

spectrum_plot = 
    D %>%    
    ggplot( ) +
    geom_point( aes(x = (mz_L+mz_R)/2, y = tot_intensity ), color = 'black', alpha = .01 ) +
    geom_rect( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = tot_intensity ), fill = 'black', alpha = .01 ) +
    geom_rect( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = tot_estimate ), fill = 'red' ) +
    geom_rangeframe() +
    theme_tufte() +
    scale_x_continuous(
        breaks = seq( min_mz, max_mz, by = 5 )
    ) 

ggsave(
    filename = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_plot.pdf',
    plot     = spectrum_plot,
    limitsize = F,
    height   = 10,
    width    = 1000 
)


# D = D %>% mutate( tot_estimate = -tot_estimate ) %>% gather( 'tot', 'value', 9:10 ) 
#     
#     
# 
# 
#     ggplot( aes(xmin = mz_L, xmax = mz_R, ymin = 0, ymax = value, fill = tot)) +
#     geom_rect() +
#     geom_rangeframe() +
#     theme_tufte()
# 
