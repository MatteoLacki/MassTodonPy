from MassTodonPy import MassTodonize
from MassTodonPy.TestScripts import substanceP, ubiquitin
from time import time
from collections import Counter, defaultdict

# mol = substanceP.copy()
mol = ubiquitin.copy()

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
                    verbose         = verbose )

def make_highchars(res, Q):
    '''Prepare the outputs of MassTodon for highcharts.'''
    precursors  = defaultdict(Counter)
    precursors_intensitites = Counter()
    fragments   = defaultdict(Counter)
    for r in res['raw_estimates']:
        for M in r['alphas']:
            if M['molType'] == 'precursor':
                ETnoDs  = M['g']
                PTRs    = Q-M['q']-M['g']
                ReactionsNo = ETnoDs + PTRs
                precursors_intensitites[M['q']] += M['estimate']
                if ReactionsNo>0:
                    precursors[ReactionsNo]['PTR'] += M['estimate']*PTRs
                    precursors[ReactionsNo]['ETnoD'] += M['estimate']*ETnoDs
            else:
                frag_type   = M['molType'][0]
                frag_number = int(M['molType'][1:])
                fragments[frag_type][frag_number] += M['estimate']
    pk = precursors.keys()
    pk.sort()
    branching_ratios = []
    for i in pk:
        PTRs = precursors[i]['PTR']
        ETnoDs = precursors[i]['ETnoD']
        try:
            branching_ratio = PTRs/float(ETnoDs)
        except ZeroDivisionError:
            branching_ratio = None
        branching_ratios.append( branching_ratio )
    probs, counts = res['basic_analysis']
    try:
        branching_ratio = probs['PTR']/probs['ETnoD']
    except ZeroDivisionError:
        branching_ratio = None
