from    MassTodonPy  import MassTodonize
import  json

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

    # todo: config should not contain an entry
    try:
        results = MassTodonize(spectrum_path= spectrum_path,
                               forPlot      = False,
                               highcharts   = True,
                               verbose      = False,
                               max_times_solve = 10,
                               **config )

        highcharts = results['highcharts']
        del results['highcharts']

        print output_path+'output.json'
        with open(output_path+'output.json', 'w') as f:
            json.dump(results,f)
    except Exception as e:
        print e
        print config

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
