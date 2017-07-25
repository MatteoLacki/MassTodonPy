library(tidyverse); library(ggthemes); library(ggiraph); library(ggthemes); library(htmlwidgets); library(rbokeh)


# path_to_data = list.files('substancesP', full.names = T)[25] 
# path_output = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/substancesP_plots' 
make_plot = function(path_to_data, path_output){
    
    experimental_settings = strsplit(path_to_data, '/')[[1]][2]
    path_output = file.path(path_output, experimental_settings)
    
    dir.create( path_output, showWarnings = F )
    
    D = read_csv( paste0(path_to_data, '/long.csv') )
    S = read_csv( paste0(path_to_data, '/short.csv') ) 
    R = read_csv( paste0(path_to_data, '/remaining_peaks.csv') ) 

    S = S %>% 
        mutate( key = paste(mz_L,mz_R)) %>%
        group_by(key) %>%
        summarize( 
            mz_L = first(mz_L),
            mz_R = first(mz_R),
            tot_intensity = round(sum(tot_intensity)),
            tot_estimate  = round(sum(tot_estimate))
        ) %>%
        mutate( mz_R = (mz_L+mz_R)/2 ) %>%
        filter( tot_intensity > 0 | tot_estimate > 0 )
    
    # Things with the same mz are aggregated: these are only artificially added 0 intensity G nodes.
    D = D %>%   
        drop_na() %>%
        mutate( key = paste(mz_L,mz_R),
                mz_L= (mz_L+mz_R)/2,
                estimate = round(estimate) ) %>%
        filter( estimate > 0 )
    
    R = R %>%
        mutate( tot_intensity = round(tot_intensity) ) %>%
        filter( tot_intensity > 0 )

    base_plot = 
        figure(
            width = 900,
            height= 500,
            title = experimental_settings
        )  %>%
        ly_rect(
            data    = D,
            xleft   = mz_L,
            xright  = mz_R,
            ybottom = 0,
            ytop    = tot_estimate_tmp,
            hover   = "m/z = [@mz_L,@mz_R] >> [@molType + @g]^{@q+} = @estimate",
            color   = 'red'
        ) %>%
        ly_rect(
            data    = S,
            xleft   = mz_L,
            xright  = mz_R,
            ybottom = 0,
            ytop    = tot_intensity,
            hover   = "m/z = [@mz_L,@mz_R]  >> @tot_intensity",
            color   = 'black'
        ) %>%
        ly_rect(
            data    = R,
            xleft   = mz_L,
            xright  = mz_R,
            ybottom = 0,
            ytop    = tot_intensity,
            hover   = "m/z = [@mz_L,@mz_R]  >> @tot_intensity",,
            color   = 'grey'
        ) %>%
        x_axis(label = "m/z") %>%
        y_axis(label = "Intensity")
    
    saveWidget(widget = base_plot, file = file.path( path_output, 'plot.html') )
    return(invisible())
}


lapply(
    list.files('substancesP', full.names = T),
    make_plot,
    path_output = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/substancesP_plots' 
)


# path_to_data = list.files('ubiquitins', full.names = T)[1]
# path_output = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiquitin_plots' 
lapply(
    list.files('ubiquitins', full.names = T),
    make_plot,
    path_output = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiquitin_plots' 
)
