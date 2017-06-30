import  os
os.environ['OMP_NUM_THREADS'] = "1"


from MassTodonPy.TestScripts.standard_datasets import substancesP, substancesP_results_macOS
from MassTodonPy import MassTodonize
import cPickle as pickle

jP      = .999
mz_prec = .05
opt_P   = .99
# cut_off = 500.0
max_times_solve = 10
multiprocesses_No = None
verbose = True
solver  = 'multiprocessing'
# solver  = 'sequential'

results = {}
for ID, mol in enumerate(substancesP):
    WH, WV = [mol['experimental_setting'][x] for x in ('WH','WV')]
    results[(ID, WH, WV)] = MassTodonize(
        fasta           = mol['fasta'],
        precursor_charge= mol['precursorCharge'],
        mz_prec         = mz_prec,
        joint_probability_of_envelope= jP,
        modifications   = mol['modifications'],
        spectrum        = mol['spectrum'],
        opt_P           = opt_P,
        # cut_off         = cut_off,
        solver          = solver,
        multiprocesses_No = multiprocesses_No,
        max_times_solve = max_times_solve,
        raw_data        = True,
        verbose         = verbose               )

# saving the data on macOS
with open('/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/substancesP_results.example', 'w') as f:
    pickle.dump(results,f)

def compare_counters(a, b, diff=.00001):
    return all( abs(a[k]-b[k]) < diff for k in set(a.keys()) | set(b.keys()))

def compare_analysis_results(A, B, diff=.00001):
    return compare_counters( A[0], B[0], diff) and compare_counters( A[1], B[1], diff)


all_the_same = {}
for key in substancesP_results_macOS:
    the_same = {}
    M = substancesP_results_macOS[key]
    N = results[key]
    for algo in ('advanced_analysis', 'basic_analysis', 'intermediate_analysis'):
        the_same[algo] = compare_analysis_results(M[algo], N[algo])
    the_same['summary']= compare_counters(M['summary'], N['summary'], 1)
    all_the_same[key] = the_same

from collections import Counter

results_counter = Counter( all_the_same[k][sk] for k in all_the_same for sk in all_the_same[k] )



if results_counter['True'] == sum(results_counter[k] for k in results_counter):
    print 'All results are within error bounds.'
    print results_counter
else:
    print 'Some results are not within error bounds.'
    print results_counter
    print [ (k,sk,all_the_same[k][sk]) for k in all_the_same for sk in all_the_same[k] if all_the_same[k][sk] == False]
    value_errors = [ results[k]['summary']['L1_error_value_error'] for k in results if results[k]['summary']['L1_error_value_error'] > 0.0 ]
    print 'L1 errors', value_errors
    print '\t Total = ', sum(value_errors)
    total_intensity = sum( results[r]['summary']['intensity_original'] for r in results)
    print 'Total Intensity =', total_intensity
    print 'Total Value Error Intensity / Total Intensity =', sum(value_errors)/float(total_intensity)
    print [results[r]['summary']['L1_error_value_error/intensity_within_tolerance'] for r in results]
    error_prone_cases = dict( (k, results[k]['raw_estimates']) for k in results if results[k]['summary']['L1_error_value_error'] > 0.0 )

    with open('errors_raw.matteo', 'w') as f:
        pickle.dump(error_prone_cases,f)
