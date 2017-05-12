from simulations import fragment, skeleton_to_table, simulate
from MassTodonPy import MassTodon
from MassTodonPy.Formulator import makeFormulas as makeF
from MassTodonPy.MatchMaker import reaction_analist_basic, reaction_analist_intermediate,  reaction_analist_upper_intermediate
import cPickle as pickle
from collections import Counter
from math import log, exp
import numpy as np

# fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'; Q = 12
fasta = 'RPKPQQFFGLM'; Q = 3

AAs = Counter(f for f in fasta)
basicAA = set(['R','H','K'])
sum([ AAs[aa] for aa in basicAA])

simulation_file = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/subP_simulation.matteo'
mol_no= 1e5; sigma = .01
L1_x = L2_x = L1_alpha = L2_alpha = 0.001
PTR = ETnoD = ETD = 1.0
pars= np.array([PTR, ETnoD, ETD])
eps = .5
I   = -log(eps)/Q**2
PTR, ETnoD, ETD = pars/sum(pars)*I


def change_key(seq, q, p, name):
    if name[0]=='c':
        name = 'c'+str(int(name[1:])-1)
    return name, q, p


def simulator_ciach(fasta, Q, mol_no, ETD,PTR, ETnoD,t=1.,return_sink=False,density=False,multiple_ETD=False):
    return simulate(('*'+fasta,Q,0,'precursor'),mol_no,ETD,PTR,ETnoD,t,return_sink,multiple_ETD,density)

molsNo = 1e4

%%time
res = simulator_ciach(fasta,Q,molsNo,ETD,PTR, ETnoD, density=False)
res = dict( ((r[3],r[1],r[2]),res[r]) for r in res)

res[('precursor',Q,0)]/molsNo


%%time
simulator_ciach(fasta,Q,1e3,ETD,PTR, ETnoD, density=True, t=0.0)


def getResults( fasta, Q, PTR, ETnoD, ETD, simulator, mol_no=1e5, sigma=.01, jP=.99, mzPrec=.05, precDigits=2, minProb=.7, cutOff=0.0, topPercent=1.0, max_times_solve = 10, L1_x=0.001, L2_x=0.001, L1_alpha=0.001, L2_alpha=0.001, verbose=False ):
    '''Perform the analysis of MassTodon on simulated data.'''
    D = simulator(fasta, Q, mol_no, intensities)
    F = dict( ( (mT, q, p), (f,bp) ) for mT,f,bp,q,p in makeF(fasta, Q, 'cz').makeMolecules(1) )
    mols    = []
    quants  = []
    for d in D:
        mType, q, p = change_key(*d)
        if mType != 'c0':
            formula, bp = F[ (mType, q, p) ]
            mols.append(   (mType, formula, bp, q, p) )
            quants.append( int(D[d]*mol_no ) )
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
    return R, FE

R,FE = getResults(fasta, Q, intsy_dict)


# dict( (r, abs(R[r][0]-R[r][1])) for r in R )
# with open(simulation_file, 'w') as f:
#     pickle.dump(sim_res, f)
