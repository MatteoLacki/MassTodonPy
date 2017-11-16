%load_ext autoreload
%autoreload 2

from MassTodonPy.Data.get_dataset import get_dataset

subP = get_dataset('substanceP')
subP.precursor.fasta
subP.precursor.q
subP.precursor.name
subP.precursor.fragmentation_type
subP.precursor
subP.instrument
subP.spectrum

data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"

substanceP

with open(data_path + "substanceP2.json", "w") as f:
    json.dump(substanceP, f)
