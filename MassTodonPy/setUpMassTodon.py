%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon   import MassTodon

fasta = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q = 8
jP= .999
mzPrec      = .05
precDigits  = 2
L2 = 0.00001
M_minProb   = .7
M = MassTodon(  fasta = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

path='/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
path_file = path+'Ubiquitin_ETD_10 ms_1071.mzXML'
cutOff = 100
topPercent = .999
M.readSpectrum(path=path_file, cutOff=cutOff, digits=precDigits, topPercent=topPercent)
M.prepare_problems(M_minProb)

%%time
res = M.run('sequential','MSE',L2=L2)
res
