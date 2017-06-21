from deconv_misc import change_key, data, sigmas2probs, probs2sigmas
from MassTodonPy import MassTodon, MassTodonize
from MassTodonPy.Formulator import make_formulas as makeF
import cPickle as pickle
from collections import Counter

sigmas = [probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629)]

fp_main = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico'
fp_in   = fp_main+'/results_Ciach/'
fp_out  = fp_main+'/results_Matteo/'

molsNo = 100000
verbose= True
solver = 'sequential'
# solver = 'multiprocessing'

sigma = sigmas[0]

with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
    res = pickle.load(f)

simulation_res = res[0]
(Q, fasta, eps, molsNo, probs), D = simulation_res

F = dict( ( (mT, q, p), (f,bp) ) for mT,f,bp,q,p in makeF(fasta, Q, 'cz').makeMolecules(1) )
mols    = []
quants  = []
for d in D:
    mType, q, p = change_key(*d)
    if mType != 'c0':
        formula, bp = F[ (mType, q, p) ]
        mols.append(   (mType, formula, bp, q, p) )
        quants.append( int(D[d]) )

RD = dict( ((mT, f, q, g), I) for (mT, f, bp, q, g), I in zip(mols, quants) )


M = MassTodon(  fasta           = fasta,
                precursor_charge= Q,
                mz_prec         = .05 )


spectrum = M.IsoCalc.makeRandomSpectrum(mols, quants, sigma, prec_digits=2)
M.read_n_preprocess_spectrum(
    spectrum = spectrum,
    opt_P    = 0.99  )
M.spectra

M.read_n_preprocess_spectrum(   spectrum  = spectrum,
                                opt_P     = .99 )


masstodon_res = MassTodonize(   fasta           = fasta,
                                precursor_charge= Q,
                                mz_prec         = .05,
                                spectrum        = spectrum,
                                opt_P           = 0.99,
                                verbose         = verbose       )

masstodon_res.keys()

masstodon_res['summary']
