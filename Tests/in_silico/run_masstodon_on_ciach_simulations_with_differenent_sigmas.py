import cPickle as pickle
from MassTodonPy import MassTodon
from MassTodonPy.Formulator import make_formulas as makeF
from collections import Counter
from math import log, exp
import numpy as np

def change_key(seq, q, p, name):
    if name[0]=='c':
        name = 'c'+str(int(name[1:])-1)
    return name, q, p


sigma=.0; jP=.99; mzPrec=.05; precDigits=2; minProb=.7; cutOff=0.0; topPercent=1.0; max_times_solve = 10; L1_x=0.001; L2_x=0.001; L1_alpha=0.001; L2_alpha=0.001; verbose=False

simulation_res

def getResults( simulation_res, sigma=.01, jP=.99, mzPrec=.05, precDigits=2, minProb=.7, cutOff=0.0, topPercent=1.0, max_times_solve = 10, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False ):

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

    M  = MassTodon( fasta           = fasta,
                    precursorCharge = Q,
                    precDigits      = precDigits,
                    jointProbability= jP,
                    mzPrec          = mzPrec )

    spectrum = M.IsoCalc.makeRandomSpectrum(mols, quants, sigma, jP, precDigits)
    M.readSpectrum(spectrum=spectrum, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
    M.prepare_problems(minProb)
    res = M.run(solver  = 'sequential',
                method  = 'MSE',
                max_times_solve = max_times_solve,
                L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                verbose = verbose )

    FE= sum(r['error'] for r in res) # fit error
    E = dict(((e['molType'], e['formula'], e['q'], e['g']), e['estimate']) for r in res for e in r['alphas'])
    RD, E = map( Counter, (RD, E))
    R = dict( (k, (RD[k], E[k])) for k in set(RD) | set(E) ) # results
    RFE = FE/sum(D[k] for k in D) # relative fit error

    ra_res = {}
    try:
        ra_res['up_int']= reaction_analist_upper_intermediate(res, Q, fasta, 0.0)
    except:
        print 'Failed upper intermediate'
    try:
        ra_res['int']   = reaction_analist_intermediate(res, Q, fasta, 0.0)
    except:
        print 'Failed intermediate'
    try:
        ra_res['base']  = reaction_analist_basic(res, Q, fasta, 0.0)
    except:
        print 'Failed basic'

    return R, FE, RFE, ra_res, res, simulation_res, spectrum

fp_in  = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Ciach/'
fp_out = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/'
# 0.0247457627119: 0.04332658,
# 0.0345762711864: 0.14815520000000004,
# molsNo=100000
# alpha=0.0247457627119
# i=10
def run(molsNo, fp_in, fp_out):
    with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
        res = pickle.load(f)
    i = 0
    sigmas = (0.0247457627119, 0.0345762711864)
    for r in res:
        for sigma in sigmas:
            print i, 'out of', len(res)*len(sigmas)
            try:
                e = getResults(simulation_res=r, sigma=sigma)
                RES = (r,e)
            except:
                RES = (r,None)
            with open( fp_out+'molsNo-'+str(molsNo)+'_sigma-'+str(sigma)+'_spec-'+str(i)+'.sadoMasto', 'wb') as handle:
                pickle.dump(RES, handle)
            i += 1
    print 'Finished with', molsNo

run(100000, fp_in, fp_out)
run(10000, fp_in, fp_out)
run(1000, fp_in, fp_out)
