import  os
os.environ['OMP_NUM_THREADS'] = "1"

from MassTodonPy.TestScripts.substancesP import substancesP
from MassTodonPy import MassTodonize
import cPickle as pickle

results = {}
for ID, mol in enumerate(substancesP):
    WH, WV = [mol['experimental_setting'][x] for x in ('WH','WV')]
    if WH == 150 and WV == 400:
        results[(ID, WH, WV)] = MassTodonize(
            fasta           = mol['fasta'],
            precursor_charge= mol['precursorCharge'],
            mz_prec         = .065,                     # mz_prec = .05
            joint_probability_of_envelope= .999,
            modifications   = mol['modifications'],
            spectrum        = mol['spectrum'],
            opt_P           = .99, # cut_off = 500.0,
            solver          = 'multiprocessing',# solver  = 'sequential'
            multiprocesses_No = None,
            max_times_solve = 10,
            raw_data        = True,
            verbose         = True )

results.keys()


results[(35, 150, 400)].keys()

results[(35, 150, 400)]['summary']
results[(35, 150, 400)]['raw_estimates'][0]['alphas']

results[(35, 150, 400)].keys()
results[(35, 150, 400)]['basic_analysis']




results[(35, 150, 400)]['raw_estimates'][0]
