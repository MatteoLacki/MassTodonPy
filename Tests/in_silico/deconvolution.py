from    deconv_misc import change_key, sigmas2probs, probs2sigmas
from    MassTodonPy import MassTodon, MassTodonize
from    MassTodonPy.Formulator import make_formulas
import  cPickle as pickle
from    collections import Counter
from    math import sqrt
import  sys
from    multiprocessing import Pool
from    itertools import repeat, product, islice
import  json

# sigmas = [probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629)]
# fp_main = sys.argv[1]
# fp_in   = fp_main+'/results_Ciach/'
# fp_out  = fp_main+'/results_Matteo/'
# verbose= False
# solver = 'sequential' # solver = 'multiprocessing'
# sigma = sigmas[0]
# molsNo = 100000
# with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
#     ciachator_res = pickle.load(f)
# simulation_res = ciachator_res[0]
# solver='sequential'
# verbose=True
def getResults( simulation_res,
                sigma,
                solver='sequential',
                verbose=False             ):
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
                                    mz_prec         = .05,
                                    spectrum        = spectrum,
                                    opt_P           = 0.999,
                                    solver          = solver,
                                    verbose         = verbose       )
        # Getting envelopes estimates
    estimates = dict(   ((e['molType'], e['formula'], e['q'], e['g']), e['estimate'])
                for r in masstodon_res['raw estimates'] for e in r['alphas']    )
    total_estimated_intensity = sum(estimates.values())
    simulated_data_cnt, estimates_cnt = map( Counter, (simulated_data_dict, estimates))
    fit_errors = dict( (k, (simulated_data_cnt[k], estimates_cnt[k]))
        for k in set(simulated_data_cnt) | set(estimates) )
        # Establishing some statistics
    stats = {}
    stats['total_fit_error_L1'] = sum( abs(real - estim) for real, estim in fit_errors.values())
    stats['total_fit_error_L2'] = sqrt(sum( (real - estim)**2 for real, estim in fit_errors.values()))
    stats['total_overestimates']= sum( max(estim - real, 0.0) for real, estim in fit_errors.values())
    stats['total_underestimates']=sum( max(real - estim, 0.0) for real, estim in fit_errors.values())
    stats['relative_to_sim_intensity_L1_error'] = stats['total_fit_error_L1']/total_simulated_intensity
    stats['relative_L1_error'] = stats['total_fit_error_L1']/(total_simulated_intensity+total_estimated_intensity)
    stats['relative_underestimates']= stats['total_underestimates']/total_simulated_intensity
    stats['relative_overestimates'] = stats['total_underestimates']/total_estimated_intensity
        # Simulated probs
    probs = {'PTR':PTR, 'ETnoD':ETnoD, 'ETD':ETD}
        # Update results
    masstodon_res['deconvolution stats']        = stats
    masstodon_res['deconvolution fit errors']   = fit_errors
    masstodon_res['probs'] = probs
    return masstodon_res

sigmas = [ probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629) ]
fp_main= sys.argv[1]
multiprocesses_No = int(sys.argv[2])

fp_in  = fp_main+'/results_Ciach'
fp_out = fp_main+'/results_Matteo'


with open(fp_main+'/data/sigmas_probs.json', 'r') as f:
    s2p = json.load(f)

sigmas2probs = dict(s2p)
probs2sigmas = dict( (b,a) for a,b in s2p )


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

        print res

        with open(fp_out+'/'+str(i), 'wb') as handle:
            pickle.dump(res, handle)
        print 'Finished with', molsNo
    except:
        OK = False
    return OK

K = len(simulated_datasets)
P = Pool(multiprocesses_No)
results = P.map(
    helper,
    zip( product(simulated_datasets, sigmas), repeat(fp_out), xrange(K) )  )
P.close()
P.join()
