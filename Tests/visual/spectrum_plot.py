from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time
from pandas import DataFrame as DF


mol = substanceP.copy()
# mol = ubiquitin.copy()

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .99,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    raw_data        = True,
                    verbose         = False )

def iBG(raw_estimates, node_type):
    for r in raw_estimates:
        SG = r['SG']
        for N in SG:
            N_D = SG.node[N]
            if N_D['type'] == node_type:
                yield N, N_D

gnodes = DF( G for _, G in iBG(res['raw_estimates'], 'G'))
inodes = DF( G for _, G in iBG(res['raw_estimates'], 'I'))
