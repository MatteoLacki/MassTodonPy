from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts.substanceP import substanceP
from MassTodonPy.TestScripts.ubiquitin  import ubiquitin
from MassTodonPy.Outputing.to_etdetective import results_to_etdetective
from time  import time

# mol = substanceP.copy()
mol = ubiquitin.copy()

res = MassTodonize( fasta           = mol['fasta'],
                    precursor_charge= mol['Q'],
                    mz_prec         = .05,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .95,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    # output_csv_path = '/Users/matteo/Documents/MassTodon/results/',
                    verbose         = False )


from collections import Counter
import numpy as np

# Q = 6
# good_data = Counter(dict( ((Q - m['q'] - m['g'], m['g']), m['estimate']) for m in some_data))
# good_data = {}
# good_data[(0,0)] = 10.
# good_data[(1,0)] = 15.
# good_data[(0,1)] = 20.
# good_data[(2,0)] = 10.
# good_data[(1,1)] = 40.
# good_data[(0,2)] = 30.
# good_data[(3,0)] = 5.
# good_data[(2,1)] = 12.
# good_data[(1,2)] = 10.
# good_data[(0,3)] = 2.
# good_data[(4,0)] = 1.
# good_data[(3,1)] = 3.
# good_data[(2,2)] = 3.
# good_data[(1,3)] = 1.
# good_data[(0,4)] = 0.
#

def is_bitonic(prev_diff, diff, tol = 0.0):
    return not (prev_diff < 0 and diff > 0)


def cantor_iter(levels):
    for a in xrange(levels):
        for i in xrange(a+1):
            yield a-i, i


def are_precursor_paths_bitonic(data, tol = 0.0):
    levels = max(a+b for a,b in data) + 1
    min_diffs   = np.zeros((levels,levels))
    intensities = np.zeros((levels,levels))
    for a,b in data:
        intensities[a,b] = data[(a,b)]
    for a,b in cantor_iter(levels):
        A = B = np.inf
        if a > 0:
            A = intensities[a,b] - intensities[a-1,b]
        if b > 0:
            B = intensities[a,b] - intensities[a,b-1]
        min_diffs[a,b] = min(A,B)
    min_diffs[0,0] = 0
    bitonic = True
    for a,b in cantor_iter(levels):
        A = B = np.inf
        if a > 0:
            diff = intensities[a,b] - intensities[a-1,b]
            prev_min_diff = min_diffs[a-1,b]
            bitonic = bitonic and is_bitonic(prev_min_diff, diff, tol)
            if not bitonic:
                print a,b, diff
                print a-1,b, prev_min_diff
                break
        if b > 0:
            diff = intensities[a,b] - intensities[a,b-1]
            prev_min_diff = min_diffs[a,b-1]
            bitonic = bitonic and is_bitonic(prev_min_diff, diff, tol)
            if not bitonic:
                print a,b, diff
                print a-1,b, prev_min_diff
                break
    return bitonic


res['raw_estimates']
Q = mol['Q']

def are_precursor_paths_bitonic_masstodon_wrapper(
        raw_estimates, Q,
        tol = 0.0, min_estimate=0.0):
    '''Wrap MassTodonPy results so that they match the precursor path checker.'''
    data = Counter(dict(((Q - alpha['q'] - alpha['g'], alpha['g']), alpha['estimate'])
           for r in res['raw_estimates']
                for alpha in r['alphas']
                    if  alpha['molType'] == 'precursor' and
                        alpha['estimate'] > min_estimate) )
    return are_precursor_paths_bitonic(data, tol = tol)


are_precursor_paths_bitonic_masstodon_wrapper( res['raw_estimates'], mol['Q'])
