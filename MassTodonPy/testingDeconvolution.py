%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import MassTodon
from    Visualization   import plot_spectrum
from    collections     import Counter
from    math            import sqrt


def simulate(   mols, sigma, fasta='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG',
                mu=1e-5, lam=0.0, M_minProb=.8, cutOff=100, topPercent=1.0, deconvMode='MSE', solverMode='sequential', jP=.999, precDigits=2, mzPrec=0.05, Q=8):
    M = MassTodon(  fasta           = fasta,
                    precursorCharge = Q,
                    precDigits      = precDigits,
                    jointProbability= jP,
                    mzPrec          = mzPrec )
    M.Forms.makeMolecules = mols.__iter__
    mol_keys, quants = zip(*mols.items())
    spectrum = M.IsoCalc.makeRandomSpectrum(mol_keys, quants, sigma, jP, precDigits)
    M.readSpectrum(spectrum=spectrum, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
    M.prepare_problems(M_minProb)
    res = M.run(solverMode, deconvMode, mu=mu, lam=lam)
    return res

mzPrec  = 0.05
x0      = 1000000
sigma   = 2*mzPrec/3
molNo   = 3
cutOff  = 100
mols = [('p','C378H629N105O118S1',76,1,0),
        ('p','C378H629N105O118S1',76,1,1),
        ('p','C378H629N105O118S1',76,1,2),
        ('p','C378H629N105O118S1',76,3,0),
        ('p','C378H629N105O118S1',76,3,1),
        ('p','C378H629N105O118S1',76,3,2) ]

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

mols = add_quants(x0, mols)


%%time
res = simulate(mols, sigma)

non_optimal = [ R for R,e,s in res if s is not 'optimal' ]
res = dict(((r['molType'], r['formula'], 76, r['q'], r['g']), r['estimate'] ) for R,e,s in res for r in R if s is 'optimal' )

error = sqrt(sum( (res[k]-mols[k])**2  for k in set(res) | set(mols)) )/x0
error
