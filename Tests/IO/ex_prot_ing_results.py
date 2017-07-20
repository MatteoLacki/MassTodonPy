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
                    raw_data        = True,
                    highcharts      = True,
                    verbose         = False )


from pandas import DataFrame as DF
import os

def i_flatten_raw_estimates(raw_estimates, threshold=0.0):
    for r in raw_estimates:
        for e in r['alphas']:
            if e['estimate'] > threshold:
                e = e.copy()
                del e['cnt']
                del e['type']
                e['status'] = r['status']
                yield e


def check_and_create_folder(path):
    if path[-1] != '/':
        path += '/'
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def write_raw_to_csv(raw_estimates, path, threshold=0.0):
    try:
        path = check_and_create_folder(path)
        res_df = DF(i_flatten_raw_estimates(raw_estimates, threshold))
        res_df = res_df[['molType','formula','q','g','estimate','status']]
        res_df = res_df.sort_values('estimate', ascending=False)
        res_df.to_csv(path_or_buf = path+'raw_estimates.csv', index = False)
        print 'Saved estimates to file.'
    except Exception as e:
        print e
        print 'Results not saved.'

# path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/IO/Output_test'
# write_raw_to_csv(res['raw_estimates'], path)
fasta = mol['fasta']


def iterate_algo(res, fasta):
    for algo in ('basic_analysis','intermediate_analysis','advanced_analysis'):
        probs, counts = res[algo]
        for k in set(probs) | set(counts):
            try:
                aa_no = int(k)
                aa = fasta[aa_no]
                L  = len(fasta)
                key = 'c'+str(k)+'-z'+str(L-k)+' bond broken on '+str(aa)
            except ValueError:
                key = k

            yield { 'Algorithm':     algo.split('_')[0],
                    'Key':          key,
                    'Probability':  probs.get(k, None),
                    'Count':        counts.get(k, None)
                   }


def write_counts_n_probs_to_csv(res, fasta, path):
    try:
        path    = check_and_create_folder(path)
        res_df  = DF(iterate_algo(res, fasta))[['Algorithm', 'Key', 'Count', 'Probability']]
        res_df.to_csv(path_or_buf = path+'counts_n_probs.csv', index = False)
        print "Saved counts'n'probs to file."
    except Exception as e:
        print e
        print 'Counts not saved.'

# write_counts_n_probs_to_csv(res, fasta, path)


def iterate_summary(summary):
    for k in summary:
        if not k in ('L1_error_value_error/original_total_intensity', 'L1_error_value_error/intensity_within_tolerance'):
            yield {'Statistic': k, 'Value': summary[k]}

def write_summary_to_csv(res, path):
    try:
        path    = check_and_create_folder(path)
        res_df  = DF(iterate_summary(res['summary']))[['Statistic', 'Value']]
        res_df.to_csv(path_or_buf = path+'summary.csv', index = False)
        print "Saved summary to file."
    except Exception as e:
        print e
        print 'Summary not saved.'

# write_summary_to_csv(res, path)
