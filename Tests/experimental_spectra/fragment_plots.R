library("jsonlite")
suppressPackageStartupMessages(library("tidyverse"))

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/substanceP_results.json'
res = read_json(path)


str(res)
