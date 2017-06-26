#TODO UPDATE THIS !!!!

from    MassTodonPy  import MassTodon
from    time import time
import  json
import  pandas as pd


def run_masstodon(  spectrum_path,
                    fasta,
                    Q,
                    modifications = {},
                    jP          =.999,
                    mzPrec      =.05,
                    precDigits  = 2,
                    M_minProb   = .7,
                    cutOff      = 100.,
                    cutOff2     = 0.0,
                    topPercent  = .999,
                    solver      = 'sequential',
                    solver_mode = 'MSE',
                    solver_max_T= 30,
                    L1_x        = 0.001,
                    L2_x        = 0.001,
                    L1_alpha    = 0.001,
                    L2_alpha    = 0.001,
                    verbose     = False ):

    '''Run MassTodon and analyze its results by basic, intermediate and upper intermediate analyzers.'''

    params = (fasta, Q, modifications, spectrum_path, jP, mzPrec, precDigits, M_minProb, cutOff, topPercent, solver_max_T, L1_x, L2_x, L1_alpha, L2_alpha)

    try:
        M = MassTodon(  fasta           = fasta,
                        precursorCharge = Q,
                        precDigits      = precDigits,
                        jointProbability= jP,
                        mzPrec          = mzPrec,
                        modifications   = modifications  )

        M.readSpectrum( path        = spectrum_path,
                        cutOff      = cutOff,
                        digits      = precDigits,
                        topPercent  = topPercent  )

        M.prepare_problems(M_minProb)
        T0_deconv = time()
        deconvolution_results = M.run(
                        solver          = solver,
                        method          = solver_mode,
                        max_times_solve = solver_max_T,
                        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
                        verbose = verbose )

        T1_deconv = time()
        T_deconv  = T1_deconv - T0_deconv

        results_analyzed = {}
        try:
            results_analyzed['base'] = M.analyze_reactions(analyzer='basic')
        except:
            print 'Basis missing'
        try:
            results_analyzed['inter'] = M.analyze_reactions(analyzer='inter')
        except:
            print 'Intermediate missing'
        try:
            results_analyzed['up_inter']= M.analyze_reactions(analyzer='up_inter')
        except:
            print 'Upper-Intermediate missing'

        if verbose:
            res = deconvolution_results, results_analyzed, params, T1_deconv
        else:
            res = deconvolution_results, results_analyzed
    except Exception as e:
        res = e, params
    return res


def get_name(key):
    return "_".join(map(str,key))

def get_subsequence(fasta, name):
    if name[0]=='p':
        return fasta
    if name[0]=='z':
        return fasta[len(fasta) - int(name[1:]):]
    else:
        return fasta[0:int(name[1:])]

def gen_data(deconvolution_results, fasta, Q):
    for r in deconvolution_results:
        for x in r['alphas']:
            name = x['molType']
            if not (name == 'precursor' and x['q'] == Q):
                f = {'seq':get_subsequence(fasta, name), 'Q':x['q'],'G':x['g'],'fragName':name, 'intensity':x['estimate'] }
                yield f

def perform_calculations(spectrum_path, output_path, file_name, config):
    '''Run MassTodonPy on a given spectrum into a given output folder.'''

    deconvolution_results, results_analyzed = \
        run_masstodon(spectrum_path = spectrum_path, **config)

    with open(output_path+'probs.json', 'w') as f:
        json.dump(results_analyzed,f)

    ParsedOutput = pd.DataFrame(gen_data(deconvolution_results, config['fasta'], config['Q']))[[2,1,0,3,4]]

    ParsedOutput.to_csv(
        path_or_buf = output_path + file_name + ".csv",
        index   = False,
        sep     = '\t' )



# spectrum_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/FRL-010513-SUBP-WH000-WV300.txt'
# fasta = 'RPKPQQFFGLM'
# Q = 3
# modifications = { 'C11':{'H':1,'O':-1,'N':1} }

# cutOff = 100; topPercent = .999; solver_max_T=30
# jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
# L1_x = L2_x = L1_alpha = L2_alpha = .001
# solver = 'sequential'; solver_mode = 'MSE'
# verbose = False
# default_config = {  'modifications' : {},
#                     'jP'            :.999,
#                     'mzPrec'        :.05,
#                     'precDigits'    : 2,
#                     'M_minProb'     : .7,
#                     'cutOff'        : 100.,
#                     'cutOff2'       : 0.0,
#                     'topPercent'    : .999,
#                     'solver'        : 'sequential',
#                     'solver_mode'   : 'MSE',
#                     'solver_max_T'  : 30,
#                     'L1_x'          : 0.001,
#                     'L2_x'          : 0.001,
#                     'L1_alpha'      : 0.001,
#                     'L2_alpha'      : 0.001,
#                     'verbose'       : False }
# run_masstodon(  spectrum_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/FRL-010513-SUBP-WH000-WV300.txt',
#                 fasta = 'RPKPQQFFGLM',
#                 Q     = 3,
#                 modifications = { 'C11':{'H':1,'O':-1,'N':1} } )
