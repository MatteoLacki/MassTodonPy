library(tidyverse)
library(jsonlite)

#################### getting parsed spectra names. these include the description of the experimental conditions
str1 = function(x) str(x, max.level=1)
load('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/info.RData')


normal_spectra_files = tools::file_path_sans_ext(list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/')[-1])
initial_run_spectra  = tools::file_path_sans_ext(list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun'))

check_which_run = function(file) ifelse( file %in% initial_run_spectra, 'initialRun', 'fullRun' )
experimentalParameters = 
  experimentalParameters %>% 
  mutate( 
    run = check_which_run(files), 
    supActive = sapply( strsplit( as.character(supActive),'\\s'), '[[',3 )
  ) 

write_json(experimentalParameters, '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubi_params.json')

#################### getting spectra
get_spectrum = function(MassTodon)
  MassTodon$spectrum    %>%
  bind_rows(.id='type') %>%
  select(mz, intensity) %>%
  na.omit()

file_names = c(list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/')[-1], list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun'))
# testing if all files are there
all(tools::file_path_sans_ext(file_names) %in% tools::file_path_sans_ext(experimentalParameters$files))

file_names = c( list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/', full.names = T)[-1],
                list.files('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/fits/initialRun', full.names = T)  )

spectra = lapply(
  file_names,
  function(x){
    load(x)
    get_spectrum(MassTodon)
  }
)

write_json(spectra, '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubi_spectra.json')


############################ analysis

# spectrum %>% plot(type='h')
# log10(spectrum$intensity) %>% hist(breaks=100)
# 
# 
# tol = .95
# i = 2
# 
# 
# show_plot = function(spectrum){
#   x = sort(spectrum$intensity, decreasing=T)
#   w = cumsum(x)/sum(spectrum$intensity) - tol > 0
#   qplot(x=log(spectrum$intensity), geom='histogram', bins = 150 ) +
#     geom_vline(xintercept = log(x[w][1]), col = 'red')
# }
# 
# plots = lapply(spectra, show_plot)
# 
# big_plot = cowplot::plot_grid(plotlist = plots, nrow=167)
# cowplot::ggsave('histograms.pdf', big_plot, limitsize = F, dpi=600, width = 100, height = 6000, units = 'mm', scale=1)


# summaries = spectra %>% lapply(summary)
# summaries[[1]]
# spectrum = spectra[[1]]


# histogram = hist(log(spectrum$intensity), breaks=100 )
#
# D = tibble(mids = histogram$mids, counts=histogram$counts) %>%
# mutate(
#   I    = exp(mids)*counts,
#   perc = I/sum(I),
#   I_cs = cumsum(I)/sum(I)
# )
#
# D %>%
# ggplot(aes(x=mids, y = counts)) +
# geom_point(aes(size=perc)) +
# geom_line() +
# theme_minimal()
#
# eps = 0.05
#
# D %>%
#   ggplot() +
#   geom_rect(
#     aes( xmin = mids-eps,
#          xmax = mids+eps,
#          ymin = 0,
#          ymax = counts)
#   )


# spectrum %>% summary
# spectrum %>% filter(intensity < 10e4) %>% plot(type='h')
#
# spectrum$intensity %>% sort(decreasing = T) %>% cumsum %>% plot(type='l')
# spectrum$intensity %>% sort(decreasing = T) %>% plot(type='l')
#
#
# library(modeest)
#
# plots = lapply(spectra, function(spectrum){
#   x = mlv(log(spectrum$intensity), method = "mfv")
#   qplot(x=log(spectrum$intensity), geom='histogram', bins=150) + geom_vline(xintercept = x$M, color = 'red')
# })
#
# spectrum
#
# hist(, plot=F)
#
# for( i in 1:length(plots)){
#   print(i)
#   plot(plots[[i]])
# }
#
#
#
# spectrum$intensity %>% sort %>% cumsum %>% plot(type='l')
# spectrum$intensity %>% sort %>% plot(type='l')
