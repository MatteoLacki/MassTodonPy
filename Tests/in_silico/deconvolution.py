from deconv_misc import change_key, data, sigmas2probs, probs2sigmas
from MassTodonPy import MassTodon, MassTodonize
from MassTodonPy.Formulator import make_formulas
import cPickle as pickle
from collections import Counter
from math import sqrt

sigmas = [probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629)]

fp_main = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico'
fp_in   = fp_main+'/results_Ciach/'
fp_out  = fp_main+'/results_Matteo/'

molsNo = 100000
# verbose= True
verbose= False
solver = 'sequential'
# solver = 'multiprocessing'

sigma = sigmas[0]

with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
    ciachator_res = pickle.load(f)


simulation_res = ciachator_res[0]


def getResults(simulation_res, verbose=False):
    '''Run MassTodon on Ciachator simulated spectra.'''

    (Q, fasta, eps, molsNo, probs), simulated_data = simulation_res

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
                                    opt_P           = 0.99,
                                    verbose         = verbose       )


        # Getting estimates
    estimates = dict(   ((e['molType'], e['formula'], e['q'], e['g']), e['estimate'])
                for r in masstodon_res['raw estimates'] for e in r['alphas']    )

    total_estimated_intensity = sum(estimates.values())

    simulated_data_cnt, estimates_cnt = map( Counter, (simulated_data_dict, estimates))

    fit_errors = dict( (k, (simulated_data_cnt[k], estimates_cnt[k]))
        for k in set(simulated_data_cnt) | set(estimates) )

    stats = {}

    stats['total_fit_error_L1'] = sum( abs(real - estim) for real, estim in fit_errors.values())
    stats['total_fit_error_L2'] = sqrt(sum( (real - estim)**2 for real, estim in fit_errors.values()))
    stats['total_overestimates']= sum( max(estim - real, 0.0) for real, estim in fit_errors.values())
    stats['total_underestimates']=sum( max(real - estim, 0.0) for real, estim in fit_errors.values())
    stats['relative_to_sim_intensity_L1_error'] = stats['total_fit_error_L1']/total_simulated_intensity
    stats['relative_L1_error'] = stats['total_fit_error_L1']/(total_simulated_intensity+total_estimated_intensity)
    stats['relative_underestimates']= stats['total_underestimates']/total_simulated_intensity
    stats['relative_overestimates'] = stats['total_underestimates']/total_estimated_intensity

masstodon_res.keys()
masstodon_res
