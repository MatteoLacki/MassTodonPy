import json
from pprint import pprint

storagePath = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/spectra.json'
with open(storagePath) as data_file:
    data = json.load(data_file)

data[0]

def extract_WH_WV(s):
    s = str(s)
    if ' ' in s:
        s = s.replace('1,5','150').replace(' ','')
    WH, WV = s.split('-')[-2:]
    WH = int(WH[2:])
    WV = int(WV[2:])
    return WH, WV

exp = data[0]
def parse_experiment(exp):
    WH, WV = extract_WH_WV(exp['instrumental_setting'][0])
    pass
