spsm = suppressPackageStartupMessages
spsm(library("jsonlite"));spsm(library("tidyverse"));spsm(library(cowplot))
len = length
analysis_results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/analyzed_results.json'

R = read_json(analysis_results_path, simpifyVector=TRUE)

parse_key = function(x) as.integer(unlist(strsplit(x,'_'))[-1])
parse_results = function(R) Map( function(x, n){
    info = parse_key(n)
    Q = info[1]; WH = info[2]; WV = info[3]
    E = lapply(x,function(y) unlist(y) %>% data_frame(name=names(.),val=., Q=Q, WH=WH, WV=WV))
    return(E) },
    R, as.list(names(R)) )

Res = parse_results(R) # Here we get ful data.
Res = lapply(1:3,function(i) bind_rows(lapply(Res,'[[',i)))
names(Res) = c('Bas', 'Int', 'UpInt')

# ________________________________________________________________________
analyze = function(D){
    DD = list(  'WH' = D %>% filter(WV==300) %>% mutate(WH=factor(WH)),
                'WV' = D %>% filter(WH==150) %>% mutate(WV=factor(WV)) )
    make_plot = function(D, par) ggplot(D) +
        geom_bar(aes_string(x=par, y='val', fill='name'), stat='identity') +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90), legend.position="top") 
    adjust_xlab = function(plot, par){
        if(par=='WV') plot = plot + xlab('Wave Velocity (Wave Height set at 1.5)')    
        if(par=='WH') plot = plot + xlab('Wave Height (Wave Velocity set at 300)')    
        return(plot)}
    shown_parameters = list(reaction        = c('reaction','no reaction'),
                        fragmentation   = c('fragmentation','no fragmentation'),
                        etnodptr        = c('ETnoD','PTR') )
    make_plots = function(par) Map(
        function(shpar) DD[[par]] %>% filter(name %in% shpar) %>% make_plot(par) %>% adjust_xlab(par),
        shown_parameters )
    return(list(WV=make_plots('WV'), WH=make_plots('WH')))
}
# Combining plots ________________________________________________________________________
plot_hahaha = function(D){
    R = analyze(D)
    list(WH=plot_grid(plotlist=R$WH, nrow=3), WV=plot_grid(plotlist=R$WV, nrow=3) )
}
plots = Map(plot_hahaha, Res)
plots$UpInt$WH