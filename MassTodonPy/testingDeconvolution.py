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
path_file = path+'dummySpec.txt'
cutOff = 100
M.readSpectrum(path_file, cutOff, precDigits)
mz, intensity = M.spectrum
len(mz)

M.prepare_problems(spectrum, M_minProb)

%%time
res = M.run('sequential','MSE',L2=L2)








readTxt(path_file)
