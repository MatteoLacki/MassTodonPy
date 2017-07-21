from MassTodonPy import MassTodonize
from MassTodonPy.Formulator import make_formulas
import cPickle as pickle
from pandas import DataFrame as DF


with open('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/data/ubiquitins.example', 'r') as h:
    ubiquitins = pickle.load(h)


# file_name = u'Ubiquitin_ETD_1428_pre_act (NCE 20)_1 ms'
file_name = u'Ubiquitin_ETD_10 ms_1071'
# The wrong spectrum.
mol =[ u for u in ubiquitins if u['experimental_setting']['files'] == file_name + '.mzXML'][0]


# This is an initial run file!!!

res = MassTodonize( fasta           = mol['fasta'],
                    # precursor_charge= mol['Q'],
                    precursor_charge= 8,
                    mz_prec         = .1,
                    joint_probability_of_envelope = .999,
                    modifications   = mol['modifications'],
                    spectrum        = mol['spectrum'],
                    opt_P           = .99,
                    solver          = 'multiprocessing',
                    multiprocesses_No = None,
                    max_times_solve = 10,
                    raw_data        = True,
                    for_plot        = True,
                    verbose         = False )


res['basic_analysis'][0]["anion_approached_cation"]
res['summary']

probs, counts = res['basic_analysis']

counts['unreacted_precursors']

[ alpha for r in res['raw_estimates'] for alpha in r['alphas']]
[ alpha for r in res['raw_estimates'] for alpha in r['alphas'] if alpha['molType'] == 'precursor']


mz, intensity = mol['spectrum']
DF({'mz':mz, 'intensity':intensity}).to_csv(
    path_or_buf = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/debugging/'+file_name +'csv',
    index = False
)

x = make_formulas('MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG', Q = 1)

list(x.makeMolecules())[0][1]
