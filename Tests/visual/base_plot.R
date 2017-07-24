library(tidyverse); library(ggthemes); library(ggiraph); library(ggthemes); library(htmlwidgets); library(rbokeh)

files = list.files('ubiquitins', full.names = T)

path_to_data = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/test/'
make_plot = function(path_to_data){
    D = read_csv( paste0(path_to_data, '/long.csv') )
    S = read_csv( paste0(path_to_data, '/short.csv') ) 
    R = read_csv( paste0(path_to_data, '/remaining_peaks.csv') ) 
    
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
    
    D %>% filter(mz_L > 1071.4, mz_R < 1071.65)
    
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
    
    saveWidget(base_plot, file = paste0('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/', path_to_data, '/plot.html'))
    return(invisible())
}

# x = lapply(files, make_plot)
# options(digits=10)
# D %>% filter(molType == 'precursor') %>% data.frame
