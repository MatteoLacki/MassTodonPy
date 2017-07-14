from MassTodonPy.TestScripts.substanceP import substanceP
from pandas import DataFrame as DF


masses, intensities = substanceP['spectrum']

DF = DF({'masses':masses, 'intensities':intensities})[['masses','intensities']]

DF.to_csv(path_or_buf='subP_spectrum.csv', index=False)
