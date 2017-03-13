%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon       import MassTodon
from    Parsers         import merge_runs
from    Visualization   import plot_spectrum
from    numpy.random    import multinomial
from    collections     import Counter
from    math            import sqrt
from    cvxopt          import matrix, spmatrix, sparse, spdiag, solvers
from    Deconvolutor    import deconvolve, Deconvolutor_Max_Flow
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

path='/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
path_file   = path+'Ubiquitin_ETD_10 ms_1071.mzXML'

cutOff      = 100
topPercent  = .999
M.readSpectrum(
    path    = path_file,
    cutOff  = cutOff,
    digits  = precDigits,
    topPercent = topPercent )

M.prepare_problems(M_minProb)



%%time
res = M.run(
    solver = 'sequential',
    method = 'MaxFlow',
    lam    = 30.,
    mu     = 0.0,
    s0_val = 2.0
)


Counter([s for a,e,s in res ])
