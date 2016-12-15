%load_ext autoreload
%autoreload
from IsotopeCalculator import isotopeCalculator, atomCnt2string
import numpy as np
from collections import Counter, defaultdict
from math import exp, fsum


fasta ='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 1
modifications = {}
ionsNo  = 10000000
P       = 2.999
isoCalc = isotopeCalculator()
digits  = 2
# molecule= {'H':1000,'C':500,'O':200,'S':20,'N':200}
molecule= {'C':100,'H':202}
masses, counts = isoCalc.randomSpectrum( molecule, ionsNo, digits, P, 0)

isoCalc.isoProbs['C']
isoCalc.isoProbs['H']

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


DF = pd.DataFrame({'mass':masses, 'intensity':counts})
DF.shape
sns.lmplot('mass', 'intensity',data=DF,fit_reg=False)
plt.show()

try:
  import cPickle as pickle
except:
  import pickle


with open('/Users/matteo/Documents/Isotoper/data/randomSpectrum.pickle', 'wb') as f:
    pickle.dump( (masses, counts), f, protocol = pickle.HIGHEST_PROTOCOL )
