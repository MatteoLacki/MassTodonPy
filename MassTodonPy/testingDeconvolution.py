%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import MassTodon
from    Visualization   import plot_spectrum
from    collections     import Counter

# from    numpy.random    import multinomial
# from    math            import sqrt
# import  numpy           as np
# import  numpy.random    as npr
# import  scipy.stats     as spy


fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8; jP = .999; mzPrec = .05; precDigits = 2
L2 = 0.00001; M_minProb = .7

M = MassTodon(  fasta = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

cutOff  = 100; topPercent = .999
x0      = 1000000
sigma   = 2*mzPrec/3

molecules = list(M.Forms.makeMolecules())
mols = molecules[0:6]
realValues = Counter()
quants = [ x0/ 2**i for i in xrange(1,len(mols)+1)]
Tq = sum(quants)
quants = [ float(q)/Tq*x0 for q in quants ]

for mol, quant in zip(mols, quants):
    _,_,_,q,g = mol
    realValues[(q,g)] = quant

spectrum = M.IsoCalc.makeRandomSpectrum(mols, quants, sigma, jP,  precDigits)


plot_spectrum( spectrum )

spectrum = M.IsoCalc.makeRandomSpectrum(mols, quants, 0.0, jP, precDigits)
plot_spectrum( spectrum )

M.readSpectrum( spectrum=(mz, counts), cutOff=100., digits=2, topPercent=1.0)
M.prepare_problems(M_minProb)

%%time
res = M.run('sequential','MSE',L2=L2)

%%time
res = M.run('sequential','MaxFlow',lam=10.)

[e for a,e,s in res]
Counter(s for a,e,s in res)

alphas, error, sol = res[0]
all(s == 'optimal' for a,e,s in res)


estimates = Counter()
for alpha in alphas:
    estimates[(alpha['q'], alpha['g'])] = alpha['estimate']

estimates
realValues

Keys = set(estimates.keys()) | set(realValues.keys())

x0
error = sum( abs( estimates[key]-realValues[key] )  for key in Keys)
error/x0 * 100
