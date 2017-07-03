import  os
old = os.environ.get('OMP_NUM_THREADS', None)
os.environ['OMP_NUM_THREADS'] = "1"

from    deconv_misc import change_key
from    MassTodonPy import MassTodon, MassTodonize
from    MassTodonPy.Formulator import make_formulas
import  cPickle as pickle
from    collections import Counter
from    math import sqrt
import  sys
from    multiprocessing import Pool
from    itertools import repeat, product, islice
import  json

# fp_main = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico'
# with open(fp_main+'/data/sigmas_probs.json', 'r') as f:
#     s2p = json.load(f)
# sigmas2probs = dict(s2p)
# probs2sigmas = dict( (b,a) for a,b in s2p )
# sigmas = [ probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629) ]
# fp_in   = fp_main+'/results_Ciach/'
# sigma = sigmas[0]
# molsNo = 100000
# with open(fp_in+'/results_molsNo-'+str(molsNo), "rb") as f:
#     ciachator_res = pickle.load(f)
# simulation_res = ciachator_res[0]
# solver='sequential'
# verbose=False
# mz_prec=.065
# opt_P=.999

def getResults( simulation_res,
                sigma,
                mz_prec     = .05,
                opt_P       = .999,
                solver      = 'sequential',
                verbose     = False     ):
    '''Run MassTodon on Ciachator simulated spectra.'''
    (Q, fasta, eps, molsNo, (PTR, ETnoD, ETD)), simulated_data = simulation_res
    formulas = dict( ( (mT, q, p), (f,bp) ) for mT,f,bp,q,p in make_formulas(fasta, Q, 'cz').makeMolecules(1) )
    total_simulated_intensity = sum(simulated_data.values())
    mols    = []
    quants  = []
    for d in simulated_data:
        mType, q, p = change_key(*d)
        if mType != 'c0':
            formula, bp = formulas[ (mType, q, p) ]
            mols.append( (mType, formula, bp, q, p) )
            quants.append( int(simulated_data[d]) )
    simulated_data_dict = dict( ((mT, f, q, g), I) for (mT, f, bp, q, g), I in zip(mols, quants) )
        # Getting spectrum
    M = MassTodon(  fasta           = fasta,
                    precursor_charge= Q,
                    mz_prec         = .05 )
    spectrum = M.IsoCalc.makeRandomSpectrum(mols, quants, sigma, prec_digits=2)
        # Running simulation
    masstodon_res = MassTodonize(   fasta           = fasta,
                                    precursor_charge= Q,
                                    mz_prec         = mz_prec,
                                    spectrum        = spectrum,
                                    opt_P           = opt_P,
                                    solver          = solver,
                                    raw_data        = True,
                                    verbose         = verbose       )
        # Getting envelopes estimates
    estimates = dict(((e['molType'], e['formula'], e['q'], e['g']), e['estimate']) for r in masstodon_res['raw_estimates'] for e in r['alphas'])
    total_estimated_intensity = sum(estimates.values())
    simulated_data_cnt, estimates_cnt = map( Counter, (simulated_data_dict, estimates))
    fit_errors = dict( (k, (simulated_data_cnt[k], estimates_cnt[k])) for k in set(simulated_data_cnt) | set(estimates) )

    stats = {}  # Establishing some statistics
    stats['total_fit_error_L1'] = sum( abs(real - estim) for real, estim in fit_errors.values())
    stats['total_fit_error_L2'] = sqrt(sum( (real - estim)**2 for real, estim in fit_errors.values()))
    stats['total_overestimates']= sum( max(estim - real, 0.0) for real, estim in fit_errors.values())
    stats['total_underestimates']=sum( max(real - estim, 0.0) for real, estim in fit_errors.values())
    stats['relative_to_sim_intensity_L1_error'] = stats['total_fit_error_L1']/total_simulated_intensity
    stats['relative_L1_error'] = stats['total_fit_error_L1']/(total_simulated_intensity+total_estimated_intensity)
    stats['relative_underestimates']= stats['total_underestimates']/total_simulated_intensity
    stats['relative_overestimates'] = stats['total_underestimates']/total_estimated_intensity
    probs = {'PTR':PTR, 'ETnoD':ETnoD, 'ETD':ETD} # Simulated probs
        # Update results
    masstodon_res['deconvolution_stats']        = stats
    masstodon_res['deconvolution_fit_errors']   = fit_errors
    masstodon_res['probs'] = probs
    del masstodon_res['raw_estimates'] # takes way too much space
    return masstodon_res



fp_main = sys.argv[1]
fp_in  = fp_main+'/results_Ciach'
fp_out = fp_main+'/results_Matteo'
multiprocesses_No = int(sys.argv[2])
with open(fp_main+'/data/sigmas_probs.json', 'r') as f:
    s2p = json.load(f)
sigmas2probs = dict(s2p)
probs2sigmas = dict( (b,a) for a,b in s2p )
sigmas = [ probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629) ]
simulated_datasets = []
for molsNo in (1000, 10000, 100000):
    with open(fp_in+'/results_molsNo-'+str(molsNo), "rb") as f:
        res = pickle.load(f)
    for r in res:
        simulated_datasets.append((r, molsNo))




def helper(helper_args):
    ((simulation_res, molsNo), sigma), fp_out, i = helper_args
    OK = True
    try:
        res = getResults(   simulation_res = simulation_res,
                            sigma = sigma,
                            solver= 'sequential',
                            verbose= False          )
        res['molsNo'] = molsNo
        res['sigma']  = sigma
        with open(fp_out+'/'+str(i)+'.matteo', 'wb') as handle:
            pickle.dump(res, handle)
        print 'Finished with', molsNo, sigma, i
    except Exception as e:
        OK = False
        print 'There is something wrong with', molsNo, sigma, i, e
    return OK


meaningful_input = product(simulated_datasets, sigmas)
all_input = zip( meaningful_input, repeat(fp_out), xrange(len(simulated_datasets)*len(sigmas)))

P = Pool(multiprocesses_No)
results = P.map( helper, all_input )
P.close()
P.join()

if old:
    os.environ['OMP_NUM_THREADS'] = old
else:
    del os.environ['OMP_NUM_THREADS']
