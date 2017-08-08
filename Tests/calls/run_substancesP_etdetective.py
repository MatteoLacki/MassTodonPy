import cPickle as pickle
import os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy                            import MassTodonize
from MassTodonPy.TestScripts.substancesP    import substancesP
from MassTodonPy.Outputing.to_etdetective   import results_to_etdetective

results = []
# mol = substancesP[0]; ID = 0
for ID, mol in enumerate(substancesP):
    fasta= mol['fasta']
    Q    = mol['precursorCharge']
    res  = MassTodonize(fasta           = fasta,
                        precursor_charge= Q,
                        mz_prec         = .05,
                        joint_probability_of_envelope = .999,
                        spectrum        = mol['spectrum'],
                        modifications   = mol['modifications'],
                        opt_P           = .99,
                        solver          = 'multiprocessing',
                        multiprocesses_No = None,
                        max_times_solve = 10,
                        raw_data        = True,
                        verbose         = False )

    mol['experimental_setting']['Q'] = Q
    WH, WV = [mol['experimental_setting'][x] for x in ('WH','WV')]
    mol['experimental_setting']['files'] = 'ID'+str(ID)+'_'+'WH'+str(WH)+'_'+'WV'+str(WV)

    etdetective_input = {
        'masstodon_output':         res,
        'experimental_setting':     mol['experimental_setting'],
        'etdetective_input':        results_to_etdetective( masstodon_results = res,
                                                            fasta = fasta,
                                                            Q = Q,
                                                            modifications = mol['modifications'] ) }

    etdetective_input['etdetective_input']['fasta'] = '*'+fasta+'*'
    results.append(etdetective_input)
    print "Dumped", ID+1, "out of", len(substancesP)

with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/Outputs/ETDetective/subsP.masstodon', 'w') as h:
    pickle.dump(results, h, protocol=pickle.HIGHEST_PROTOCOL)

# from collections import Counter
# [s['experimental_setting'] for s in substancesP]
