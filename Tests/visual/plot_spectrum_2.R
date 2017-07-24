library(tidyverse); library(ggthemes); library(ggiraph); library(ggthemes); library(htmlwidgets); library(rbokeh)

D = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_long.csv')
S = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_short.csv') 
R = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_remaining_peaks.csv') 

S = S %>% 
    mutate( key = paste(mz_L,mz_R)) %>%
    group_by(key) %>%
    summarize( 
        mz_L = first(mz_L),
        mz_R = first(mz_R),
        tot_intensity = sum(tot_intensity),
        tot_estimate  = sum(tot_estimate)
    ) %>%
    mutate(
        mz_R = (mz_L+mz_R)/2
    )

# Things with the same mz are aggregated: these are only artificially added 0 intensity G nodes.
D = D %>%   
    drop_na() %>%
    mutate( 
        key = paste(mz_L,mz_R),
        mz_L= (mz_L+mz_R)/2
    )  
    
# figure(
#     width = 900,
#     height= 500
# )  %>%
# ly_rect(
#     data    = S,
#     xleft   = mz_L,
#     xright  = mz_R,
#     ybottom = 0,
#     ytop    = tot_intensity,
#     alpha   = .5,
#     hover   = "Total Intensity: @tot_intensity"
# )

base_plot = 
    figure(
        width = 900,
        height= 500
    )  %>%
    ly_rect(
        data    = D,
        xleft   = mz_L,
        xright  = mz_R,
        ybottom = 0,
        ytop    = tot_estimate_tmp,
        hover   = "q = @q, g = @g, @molType, @estimate",
        color   = 'red'
    ) %>%
    ly_rect(
        data    = S,
        xleft   = mz_L,
        xright  = mz_R,
        ybottom = 0,
        ytop    = tot_intensity,
        hover   = "Total Intensity: @tot_intensity",
        color   = 'black'
    ) %>%
    ly_rect(
        data    = R,
        xleft   = mz_L,
        xright  = mz_R,
        ybottom = 0,
        ytop    = tot_intensity,
        hover   = "Total Intensity: @tot_intensity",
        color   = 'grey'
    )   
        

saveWidget(base_plot, file = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/bokeh_test.html')


# # bokeh_base = 
#     figure(
#         width = 900,
#         height= 500
#     )  %>%
#     ly_rect(
#         xleft   = mz_L,
#         xright  = mz_R,
#         ybottom = 0,
#         ytop    = tot_intensity_tmp,
#         alpha   = .5,
#         hover   = "Total Intensity: @tot_intensity_tmp",
#         data    = S
#     )
# 
# %>%
#     ly_rect(
#         xleft   = (mz_R+mz_L)/2,
#         xright  = mz_R,
#         ybottom = 0,
#         ytop    = tot_estimate_tmp,
#         color   = 'red',
#         alpha   = .5,
#         hover   = "q = @q, g = @g, @molType, @estimate" 
#     )
# 
# 
# 
# D %>% filter(key == '1049.55 1049.63') %>% data.frame
# options(scipen=999)
# 
# 
# 
# 
# 
# 
# 
# 
# D %>%   
#     drop_na() %>%
#     mutate( mz = (mz_L+mz_R)/2) %>%
#     group_by( mz ) %>%
#     summarize( 
#         test_total_estimate_tmp = sum(estimate),
#         total_estimate_tmp      = first(tot_estimate_tmp)
#     ) %>% 
#     data.frame
#     
# bokeh_base = 
#     D %>%    
#     drop_na() %>%
#     mutate( estimate = round(estimate) ) %>%
#     filter( mz_L > 1060, mz_R < 1080 ) %>%
#     figure(
#         width = 900,
#         height= 500
#     )  %>%
#     ly_rect(
#         xleft   = mz_L,
#         xright  = (mz_R+mz_L)/2,
#         ybottom = 0,
#         ytop    = tot_intensity_tmp,
#         alpha   = .5,
#         hover   = "Total Intensity: @tot_intensity_tmp" 
#     ) %>%
#     ly_rect(
#         xleft   = (mz_R+mz_L)/2,
#         xright  = mz_R,
#         ybottom = 0,
#         ytop    = tot_estimate_tmp,
#         color   = 'red',
#         alpha   = .5,
#         hover   = "q = @q, g = @g, @molType, @estimate" 
#     )
# 
# 
# 

# 
# 
