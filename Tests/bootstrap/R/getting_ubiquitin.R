library(tidyverse)
library(jsonlite)

#################### getting parsed spectra names. these include the description of the experimental conditions
str1 = function(x) str(x, max.level=1)
load('/Users/matteo/Dropbox/MassTodon/ProcessedData/ubiquitin/Orbi_2014_Dec/info.RData')
write_json(as.list(experimentalParameters), '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubi_params.json')

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

