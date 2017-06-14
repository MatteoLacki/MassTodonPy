from    MassTodonPy     import MassTodonize
from    numpy.random    import multinomial
import  pandas      as pd
import  cPickle     as pickle
import  json
from    bootstrap_misc import MassTodon_bootstrap



def process_ra(info, ra_type_row, prob_count):
    row_ra = {}
    for k in info:
        k_tmp = ra_type_row + "_" + prob_count + "_" + str(k).replace(" ","_")
        row_ra[k_tmp] = info[k]
    return row_ra


def get_row(res):
    row = {}
    ra_type_row = { 'basic analysis':               'basic',
                    'intermediate analysis':        'inter',
                    'upper intermediate analysis':  'up_inter' }
    for ra_type in ra_type_row:
        for i, what in enumerate(('prob','count')):
            row.update( process_ra( res[ra_type][i], ra_type_row[ra_type], what ) )
    row.update(res['summary'])
    return row


# experiments = substancesP
def MassTodon_real(experiments):
    '''Run MassTodon on real experiments.'''
    for ID, exp in enumerate(experiments):
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
        row = get_row(Results)
        row['WH'] = exp['experimental_setting']['WH']
        row['WV'] = exp['experimental_setting']['WV']
        row['ID'] = ID
        yield row


data_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/substanceP_spectra_parsed.cPickle'
with open(data_path, 'r') as f:
    substancesP = pickle.load(f)

results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/RESULTS_REAL/real.csv'

DF = pd.DataFrame(MassTodon_real(substancesP))
DF.to_csv(path_or_buf = results_path, index = False)
