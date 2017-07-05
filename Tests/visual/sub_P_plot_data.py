from    MassTodonPy import MassTodonize
import  pandas      as pd
import  cPickle     as pickle
import  json

sub_P_parsed = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(sub_P_parsed, 'r') as f:
    experiments = pickle.load(f)

path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/visual/data/'
experimental_settings = []

for i, exp in enumerate(experiments):
    Results = MassTodonize(
        fasta           = exp['fasta'],
        precursor_charge= exp['precursorCharge'],
        joint_probability_of_envelope = .990,
        mz_prec         = .05,
        modifications   = exp['modifications'],
        spectrum        = exp['spectrum'],
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

    WH = exp['experimental_setting']['WH']
    WV = exp['experimental_setting']['WV']

    experimental_settings.append((WH,WV))

    pd.DataFrame(Results['short data to plot']).to_csv(
        path_or_buf = path+str(i)+'_shortData.csv',
        index       = False)

    pd.DataFrame(Results['long data to plot']).to_csv(
        path_or_buf = path+str(i)+'_longData.csv',
        index       = False)


with open(path+'experimental_settings.json', 'w') as f:
    json.dump(experimental_settings,f)
