from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from time  import time

mol = substanceP.copy()
# mol = ubiquitin.copy()

jP      = .999
mzPrec  = .05
opt_P   = .99
max_times_solve = 10
multiprocesses_No = None
verbose = True
solver  = 'multiprocessing'
# solver  = 'sequential'

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = mzPrec,
                    joint_probability_of_envelope = jP,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = opt_P,
                    solver          = solver,
                    multiprocesses_No = multiprocesses_No,
                    max_times_solve = max_times_solve,
                    raw_data        = False,
                    highcharts      = True,
                    verbose         = verbose )


res.keys()



from pandas import DataFrame as DF

def i_flatten_raw_estimates(raw_estimates):
    for r in raw_estimates:
        for e in r['alphas']:
            e = e.copy()
            del e['cnt']
            del e['type']
            e['status'] = r['status']
            yield e


def write_raw_to_csv(raw_estimates, path):
    try:
        res_df = DF(i_flatten_raw_estimates(raw_estimates))
        res_df = res_df[['molType','formula','q','g','estimate','status']]
        res_df = res_df.sort_values('estimate', ascending=False)
        res_df.to_csv(path_or_buf = path, index = False)
        print 'Saved estimates to file.'
    except:
        print 'Results not saved.'
