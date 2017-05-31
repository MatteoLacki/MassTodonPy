from    MassTodonPy     import MassTodonize
from    numpy.random    import multinomial

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


# exp = experiments[0]
def MassTodon_bootstrap(exp, ions_no, bootstrap_size, ID, WH, WV):
    '''Perform bootstrap analysis of MassTodon results.'''
    mzs, intensities = exp['spectrum']
    # sim_intensities = multinomial(ions_no, intensities/intensities.sum())
    i = 0
    for sim_intensities in multinomial(ions_no, intensities/intensities.sum(), bootstrap_size):
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
        row = get_row(Results)
        row['WH'], row['WV'], row['ID'] = WH, WV, ID
        yield row
        if divmod(i, 10)[1] == 0:
            print 'Finished',i+1,'out of',bootstrap_size,'.'
        i += 1
