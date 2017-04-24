import  json
import  numpy as np
from    time import time
from    MassTodonPy  import MassTodon
import  cPickle as pickle
from    MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate, reaction_analist_upper_intermediate
from    collections import Counter

result_path_specific = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/parsed_sub_P.matteo'

with open(result_path_specific,'r') as fp:
    experiments = pickle.load(fp)

experiments
