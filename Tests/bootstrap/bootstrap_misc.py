from    MassTodonPy     import MassTodonize
from    numpy.random    import multinomial
import  numpy as np

def process_ra(info, ra_type_row, prob_count):
    row_ra = {}
    for k in info:
        k_tmp = ra_type_row + "_" + prob_count + "_" + str(k).replace(" ","_")
        row_ra[k_tmp] = info[k]
    return row_ra

# Results['summary']
# Results.keys()
# res = Results
def get_row(res):
    row = {}
    ra_type_row = { 'basic analysis':               'basic',
                    'intermediate analysis':        'inter',
                    'advanced analysis':            'up_inter' }
    for ra_type in ra_type_row:
        for i, what in enumerate(('prob','count')):
            row.update( process_ra( res[ra_type][i], ra_type_row[ra_type], what ) )
    row.update(res['summary'])
    return row


# exp = experiments[0]
def MassTodon_bootstrap(exp, ions_no, bootstrap_size, ID, WH, WV):
    '''Perform bootstrap analysis of MassTodon results.'''
    mzs, intensities = exp['spectrum']
    # sim_intensities = multinomial(ions_no, intensities/intensities.sum()).astype(np.float) * intensities.sum()/float(ions_no)
    i = 0
    for sim_intensities in multinomial(ions_no, intensities/intensities.sum(), bootstrap_size).astype(np.float) * intensities.sum()/float(ions_no):
        Results = MassTodonize(
            fasta           = exp['fasta'],
            precursor_charge= exp['precursorCharge'],
            mz_prec         = .05,
            opt_P           = .999,
            modifications   = exp['modifications'],
            spectrum        = (mzs, sim_intensities),
            solver          = 'sequential',
            max_times_solve = max_times_solve,
            raw_data        = False,
            verbose         = False )
        row = get_row(Results)
        row['WH'], row['WV'], row['ID'] = WH, WV, ID
        yield row
        if divmod(i, 10)[1] == 0:
            print 'Finished',i+1,'out of',bootstrap_size,'.'
        i += 1


for sim_intensities in multinomial(ions_no, intensities/intensities.sum(), bootstrap_size).astype(np.float) * intensities.sum()/float(ions_no):
