library(tidyverse); library(ggthemes); library(ggiraph); library(ggthemes); library(htmlwidgets); library(rbokeh)

D = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_long.csv')
S = read_csv('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubi_data_short.csv') 

S = S %>% 
    mutate( key = paste(mz_L,mz_R)) %>%
    group_by(key) %>%
    summarize( 
        mz_L = first(mz_L),
        mz_R = first(mz_R),
        tot_intensity = sum(tot_intensity),
        tot_estimate  = sum(tot_estimate)
    ) # Things with the same mz are aggregated: these are only artificially added 0 intensity G nodes.


D %>%   
drop_na() %>%
mutate( key = paste(mz_L,mz_R)) %>%
group_by(key)



options(scipen=999)



B = left_join( S, D, by='mz' ) 

D %>% filter(mz == 272.185) %>% data.frame



D %>%   
    drop_na() %>%
    mutate( mz = (mz_L+mz_R)/2) %>%
    group_by( mz ) %>%
    summarize( 
        test_total_estimate_tmp = sum(estimate),
        total_estimate_tmp      = first(tot_estimate_tmp)
    ) %>% 
    data.frame
    
bokeh_base = 
    D %>%    
    drop_na() %>%
    mutate( estimate = round(estimate) ) %>%
    filter( mz_L > 1060, mz_R < 1080 ) %>%
    figure(
        width = 900,
        height= 500
    )  %>%
    ly_rect(
        xleft   = mz_L,
        xright  = (mz_R+mz_L)/2,
        ybottom = 0,
        ytop    = tot_intensity_tmp,
        alpha   = .5,
        hover   = "Total Intensity: @tot_intensity_tmp" 
    ) %>%
    ly_rect(
        xleft   = (mz_R+mz_L)/2,
        xright  = mz_R,
        ybottom = 0,
        ytop    = tot_estimate_tmp,
        color   = 'red',
        alpha   = .5,
        hover   = "q = @q, g = @g, @molType, @estimate" 
    )



saveWidget(bokeh_base, file = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/bokeh_test.html')    


