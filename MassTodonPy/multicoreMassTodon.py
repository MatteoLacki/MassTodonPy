%load_ext autoreload
%load_ext line_profiler
%autoreload
from    MassTodon   import MassTodon
import  numpy       as np
from    collections import Counter, defaultdict
from    Parsers     import ParseMzXML, trim_spectrum
from    cvxopt      import matrix, spmatrix, sparse, spdiag, solvers
import  cPickle     as pickle
from    Visualization       import plot_spectrum, plot_deconvolution_graph
import  matplotlib.pyplot   as plt

spectrum_intensity_cut_off = 1000
path = '/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/data/'
spectrum= ParseMzXML(   path+'Ubiquitin_ETD_10 ms_1071.mzXML',
                        cut_off_intensity = spectrum_intensity_cut_off )
spectrum= trim_spectrum(spectrum) # plot_spectrum(spectrum, 0, 4000)
fasta   = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
Q   = 8
jP  = .999
modifications = {}
mzPrec      = .05
precDigits  = 2
L2 = 0.00001
M_minProb   = .7
M = MassTodon(  fasta = fasta,
                precursorCharge = Q,
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec )

M.prepare_problems(spectrum, M_minProb)

%%time
res = M.run('sequential','MSE',L2=L2)

alphas, error, status = res[0]
Counter(s for a,e,s in res)

##### in silico spectrum
N  = 100000
jP = .9999
masses, cnts = M.IsoCalc.randomSpectrum(ionsNo=N, atomCnt_str='C100H200', jointProb=jP)
masses, cnts

M.IsoCalc.addNoise(masses, cnts)
res[0][0][0]

res2 = defaultdict(list)

res
for r, e, s in res:
    print r['q'], r['molType']
    # res2[(r['q'], r['molType'])].append( (r,e) )

list(M.Forms.makeMolecules())
