from    MassTodon       import MassTodon
from    collections     import Counter
from    math            import sqrt
import numpy as np
import scipy.stats as ss

def analyze_res(res, mols):
    x0 = sum(mols[k] for k in mols)
    non_optimal = [ R for R,e,s in res if s is not 'optimal' ]
    res   = dict(((r['molType'], r['formula'], 76, r['q'], r['g']), r['estimate'] ) for R,e,s in res for r in R if s is 'optimal' )
    error = sqrt(sum( (res[k]-mols[k])**2  for k in set(res) | set(mols)) )/x0
    return error, non_optimal


def simulate(   mols, sigma, fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
                mu=1e-5, lam=0.0, M_minProb=.8, cutOff=100, topPercent=1.0, deconvMode='MSE', solverMode='sequential', jP=.999, precDigits=2, mzPrec=0.05, Q=8):
    M = MassTodon(  fasta           = fasta,
                    precursorCharge = Q,
                    precDigits      = precDigits,
                    jointProbability= jP,
                    mzPrec          = mzPrec )
    M.Forms.makeMolecules = mols.__iter__ # this is a witty way to avoid checking all other possible molecules
    mol_keys, quants = zip(*mols.items())
    spectrum = M.IsoCalc.makeRandomSpectrum(mol_keys, quants, sigma, jP, precDigits)
    M.readSpectrum(spectrum=spectrum, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
    M.prepare_problems(M_minProb)
    res = M.run(solverMode, deconvMode, mu=mu, lam=lam)
    return analyze_res(res,mols)


def finiteGeometric(x0, molNo):
    quants  = [ x0/ 2**i for i in xrange(molNo)]
    Tq      = sum(quants)
    quants  = [ float(q)/Tq*x0 for q in quants ]
    return quants


def add_quants(x0, mols):
    quants = finiteGeometric(x0, len(mols))
    realValues = Counter()
    for mol, quant in zip(mols,quants):
        realValues[mol]=quant
    return realValues


def getMols(molsNo):
    return [('p','C378H629N105O118S1',76,1,i) for i in xrange(molsNo)]


def get_p_gauss(sigma, mzPrec=0.05):
    N = ss.norm(loc=0, scale=sigma)
    return N.cdf(mzPrec)-N.cdf(-mzPrec)


def list_for_shell_script(arg):
    return " ".join( str(x) for x in arg)


## These are the standard deviations we test.
# sigmas = np.linspace(0.01,.3,60)
# array([ 0.01      ,  0.01491525,  0.01983051,  0.02474576,  0.02966102,
#         0.03457627,  0.03949153,  0.04440678,  0.04932203,  0.05423729,
#         0.05915254,  0.0640678 ,  0.06898305,  0.07389831,  0.07881356,
#         0.08372881,  0.08864407,  0.09355932,  0.09847458,  0.10338983,
#         0.10830508,  0.11322034,  0.11813559,  0.12305085,  0.1279661 ,
#         0.13288136,  0.13779661,  0.14271186,  0.14762712,  0.15254237,
#         0.15745763,  0.16237288,  0.16728814,  0.17220339,  0.17711864,
#         0.1820339 ,  0.18694915,  0.19186441,  0.19677966,  0.20169492,
#         0.20661017,  0.21152542,  0.21644068,  0.22135593,  0.22627119,
#         0.23118644,  0.23610169,  0.24101695,  0.2459322 ,  0.25084746,
#         0.25576271,  0.26067797,  0.26559322,  0.27050847,  0.27542373,
#         0.28033898,  0.28525424,  0.29016949,  0.29508475,  0.3       ])
# list_for_shell_script(sigmas)

## These are the corresponding probabilities of them falling into the tollerance interval [ -0.05, 0.05 ]
# get_p_gauss(sigmas, mzPrec=0.05)
# array([ 0.99999943,  0.99919849,  0.98831003,  0.95667342,  0.90814916,
#         0.8518448 ,  0.79452072,  0.73981509,  0.68929589,  0.64340622,
#         0.60204112,  0.56485823,  0.53143491,  0.50134371,  0.4741858 ,
#         0.44960344,  0.4272823 ,  0.40694922,  0.388368  ,  0.37133488,
#         0.35567401,  0.34123345,  0.32788163,  0.31550433,  0.30400205,
#         0.29328792,  0.28328577,  0.27392868,  0.2651576 ,  0.25692031,
#         0.24917047,  0.24186683,  0.2349726 ,  0.22845482,  0.22228395,
#         0.21643338,  0.21087914,  0.20559956,  0.20057498,  0.19578757,
#         0.19122111,  0.18686079,  0.18269311,  0.17870568,  0.17488715,
#         0.17122712,  0.16771598,  0.16434491,  0.16110575,  0.15799095,
#         0.15499353,  0.15210704,  0.14932545,  0.14664319,  0.14405505,
#         0.14155621,  0.13914213,  0.1368086 ,  0.13455168,  0.13236767])
