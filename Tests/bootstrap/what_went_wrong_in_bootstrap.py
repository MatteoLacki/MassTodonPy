import  pandas      as pd
import  cPickle     as pickle
import  json
from    bootstrap_misc  import  MassTodon_bootstrap
from    numpy.random    import  multinomial
import  numpy           as      np
from    collections     import  Counter
from    bootstrap_misc  import  get_row, process_ra
from    MassTodonPy     import  MassTodonize

data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

ions_no = 10**5
trouble_cases = filter(lambda x: x['experimental_setting']['WH']==150 and x['experimental_setting']['WV']==300,substancesP)
exp = trouble_cases[0]

WH = exp['experimental_setting']['WH']
WV = exp['experimental_setting']['WV']

mzs, intensities = exp['spectrum']
sim_intensities  = multinomial(ions_no, intensities/intensities.sum())
sim_intensities  = sim_intensities.astype(np.float)

Results = MassTodonize(
    fasta           = exp['fasta'],
    precursor_charge= exp['precursorCharge'],
    joint_probability_of_envelope = .990,
    mz_prec         = .05,
    modifications   = exp['modifications'],
    spectrum        = (mzs, sim_intensities),
    opt_P           = .99,
    min_prob_of_envelope_in_picking = .7,
    solver          = 'sequential',
    method          = 'MSE',
    max_times_solve = 10,
    L1_x            = .001,
    L2_x            = .001,
    L1_alpha        = .001,
    L2_alpha        = .001,
    verbose         = False )

raw = Results['raw estimates']
[ k for r in raw for k in r['alphas'] ]

SG = raw[0]['SG']
SG.nodes(data=True)
SG.edges(data=True)


row = get_row(Results)
row['WH'], row['WV'], row['ID'] = WH, WV, ID
yield row
if divmod(i, 10)[1] == 0:
    print 'Finished',i+1,'out of',bootstrap_size,'.'
i += 1


DF = pd.DataFrame( MassTodon_bootstrap(exp, ions_no, 1, 1, WH, WV) )

DF

final_path = results_path+'WH-'+str(WH)+'_WV-'+str(WV)+'_ID-'+str(ID)+'.csv'
DF.to_csv(path_or_buf = final_path, index = False)
print
print 'Dumped',ID,'out of',len(substancesP),'.'
print

#TODO: Lesson learned: keep only one naming convention. Always.
