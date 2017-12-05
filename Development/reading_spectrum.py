%load_ext autoreload
%autoreload 2

import csv

from MassTodonPy.Data.get_dataset import get_dataset
from MassTodonPy.Spectra.Read import read_spectrum

path = '/Users/matteo/Downloads/spektra_Posen/H4-qb.mzML'
spectrum = read_spectrum(path)

# sum of intensities
spectrum.total_intensity()

# full iteration
for mz, intensity in spectrum:
    print(mz, intensity)

# subsetting
for mz, intensity in spectrum[600, 800]:
    print(mz, intensity)

# visualization

# spectrum.plot()  # this is how it is supposed to look

# So we go to R...
path_to_csv = '/Users/matteo/Documents/MassTodon/MassTodonPy/Development/Bokeh/H4-qb.csv'
with open(path_to_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    for mz, intensity in spectrum:
        writer.writerow([mz, intensity])
