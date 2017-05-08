import cPickle as pickle
from MassTodonPy import MassTodon
from MassTodonPy.Formulator import makeFormulas as makeF
from MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate,  reaction_analist_upper_intermediate
from collections import Counter
from math import log, exp
import numpy as np


def change_key(seq, q, p, name):
    if name[0]=='c':
        name = 'c'+str(int(name[1:])-1)
    return name, q, p


# sigma=.01; jP=.99; mzPrec=.05; precDigits=2; minProb=.7; cutOff=0.0; topPercent=1.0; max_times_solve = 10; L1_x=0.001; L2_x=0.001; L1_alpha=0.001; L2_alpha=0.001; verbose=False

def getResults( simulation_res, sigma=.01, jP=.99, mzPrec=.05, precDigits=2, minProb=.7, cutOff=0.0, topPercent=1.0, max_times_solve = 10, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False ):
    '''I dunwanna write no docs!'''
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
        ra_res['up_int']= reaction_analist_upper_intermediate(res, Q, fasta)
    except:
        print 'Failed upper intermediate'
    try:
        ra_res['int']   = reaction_analist_intermediate(res, Q, fasta)
    except:
        print 'Failed intermediate'
    try:
        ra_res['base']  = reaction_analist_basic(res, Q, fasta)
    except:
        print 'Failed basic'    
    return R, FE, RFE, ra_res, res, simulation_res, spectrum


fp_in  = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Ciach/'
fp_out = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/'

for molsNo in (100000, 10000, 1000):
    with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
        res = pickle.load(f)

    MassTodonRes = [ getResults(r) for r in res ]

    with open(fp_out+'results_molsNo-'+str(molsNo), "w") as f:
        res = pickle.load(f)
