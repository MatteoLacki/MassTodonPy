%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import MassTodon
from    Parsers         import merge_runs
from    Visualization   import plot_spectrum
from    numpy.random    import multinomial
from    collections     import Counter
from    math            import sqrt
import  numpy as np

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8; jP = .999; mzPrec = .05; precDigits = 2
L2 = 0.00001; M_minProb = .7

M = MassTodon(  fasta = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )
cutOff = 100; topPercent = .999
molecules = list(M.Forms.makeMolecules())

ave_mz, ave_intensity = M.IsoCalc.isoEnvelope( molecules[0][1], jP=.999, q=0, g=0, precDigits=2 )
mols = molecules[0:6]


def get_envelopes(mols, quants):
    for mol, quant in zip(mols, quants):
        _, atomCnt_str, _, q, g = mol
        ave_mz, ave_intensity = M.IsoCalc.isoEnvelope(atomCnt_str=atomCnt_str, jP=.999, q=q, g=g, precDigits=2)
        ave_intensity = quant * ave_intensity
        yield ave_mz, ave_intensity

x0 = 1000000
quants = [ x0/ 2**i for i in xrange(1,len(mols)+1)]
mz, intensity   = reduce(merge_runs, get_envelopes(mols, quants))

realValues = Counter()
for mol, quant in zip(mols, quants):
    _,_,_,q,g = mol
    realValues[(q,g)] = quant

probs   = intensity/sum(intensity)
counts  = np.array( multinomial(x0, probs), dtype='f')
# plot_spectrum( (mz, counts) )

M.readSpectrum( spectrum=(mz, counts), cutOff=100., digits=2, topPercent=1.0)
M.prepare_problems(M_minProb)

%%time
res = M.run('sequential','MSE',L2=L2)

len(res)
alphas, error, sol = res[0]
all(s == 'optimal' for a,e,s in res)

estimates = Counter()
for alpha in alphas:
    estimates[(alpha['q'], alpha['g'])] = alpha['estimate']

Keys = set(estimates.keys()) | set(realValues.keys())

x0
error = sum( abs( estimates[key]-realValues[key] )  for key in Keys)
error/x0 * 100
