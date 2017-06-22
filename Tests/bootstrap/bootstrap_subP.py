import  pandas      as pd
import  cPickle     as pickle
import  json
from    bootstrap_misc import MassTodon_bootstrap

data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

ions_no         = 10**6
bootstrap_size  = 1
results_path    = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_CSV/'

ID, exp = 0, substancesP[0]

# for ID, exp in enumerate(substancesP):
WH = exp['experimental_setting']['WH']
WV = exp['experimental_setting']['WV']
DF = pd.DataFrame( MassTodon_bootstrap(exp, ions_no, bootstrap_size, ID, WH, WV) )
final_path = results_path+'WH-'+str(WH)+'_WV-'+str(WV)+'_ID-'+str(ID)+'.csv'
DF.to_csv(path_or_buf = final_path, index = False)
print
print 'Dumped',ID,'out of',len(substancesP),'.'
print

#TODO: Lesson learned: keep only one naming convention. Always.
