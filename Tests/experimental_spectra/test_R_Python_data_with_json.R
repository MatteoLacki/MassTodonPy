suppressPackageStartupMessages(library("jsonlite"))

DATA = list(a='dupa','b'=10,c=rnorm(10))
DATA

toJSON(DATA, digits=10)
write_json(DATA, path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/test_json_of_a_bitch.json', digits=10)
