import cPickle as pickle
from pandas import DataFrame as DF

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/Data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)

csv_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/ubiqutins/'

for ubi in ubiquitins:
    mass, intensity = ubi['spectrum']
    file_name = ubi['experimental_setting']['files'].split('.')[0]
    DF({'mass': mass, 'intensity':intensity }).to_csv(
        path_or_buf= csv_path + file_name + '.csv', index=False )
