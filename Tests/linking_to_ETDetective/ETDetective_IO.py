import cPickle as pickle

path_to_files = '/Users/matteo/Downloads/sample/'

with open(path_to_files+'example_data', 'r') as f:
    example_data = pickle.load(f)

with open(path_to_files+'example_precursor', 'r') as f:
    example_precursor = pickle.load(f)


from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter

mol = substanceP.copy()
cutOff = 100; topPercent = .999; max_times_solve=30
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursorCharge = mol['Q'],
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = mol['modifications']  )
M.readSpectrum( spectrum        = mol['spectrum'],
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent  )
M.prepare_problems(M_minProb)
Results = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )


def get_subsequence(fasta, name, modifications):
    if modifications == {'C11': {'H': 1, 'N': 1, 'O': -1}}:
        suffix = '*'
    else:
        suffix = ''
    if name[0]=='p':
        return '*'+fasta+suffix
    if name[0]=='z':
        return fasta[len(fasta) - int(name[1:]):]+suffix
    else:
        return '*'+fasta[0:int(name[1:])]


def flatten_results(results, fasta, modifications):
    for subproblem in results:
        for r in subproblem['alphas']:
            subseq = get_subsequence(fasta, r['molType'], modifications)
            yield (subseq, r['q'], r['g'], r['molType']), r['estimate']


def gen_ETDetective_inputs(fasta, Q, results, modifications={}):
    '''Get inputs for ETDetective.'''
    data = {}
    for r in flatten_results(results, fasta, modifications):
        (subseq, q, g, molType), estimate = r
        if q == Q and g == 0 and molType=='precursor':
            precursor = subseq, q, g, molType
        else:
            data[(subseq, q, g, molType)] = estimate
    return precursor, data

gen_ETDetective_inputs(mol['fasta'], mol['Q'], Results, mol['modifications'])

example_data
example_precursor
