from MassTodonPy import MassTodon
from MassTodonPy.TestScripts import substanceP, ubiquitin
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter


mol = substanceP.copy()

cutOff = 100; topPercent = .999; max_times_solve=30
jP=.999; mzPrec=.05; precDigits=2; M_minProb=.7
L1_x = L2_x = L1_alpha = L2_alpha = .001
solver = 'sequential'; method  = 'MSE'
verbose = False

M = MassTodon(  fasta           = mol['fasta'],
                precursorCharge = mol['Q'],
                precDigits      = precDigits,
                jointProbability= jP,
                mzPrec          = mzPrec,
                modifications   = mol['modifications']  )

M.readSpectrum( spectrum        = mol['spectrum'],
                cutOff          = cutOff,
                digits          = precDigits,
                topPercent      = topPercent  )

M.run(  solver  = 'sequential',
        method  = 'MSE',
        min_prob_per_molecule = M_minProb,
        max_times_solve = max_times_solve,
        L1_x=L1_x, L2_x=L2_x, L1_alpha=L1_alpha, L2_alpha=L2_alpha,
        verbose = verbose )

M.summarize_results()

A, B = M.gen_ETDetective_inputs()
sum(B[k] for k in B)

sum( k['estimate'] for r in M.res for k in r['alphas'] )


M.res[0]['small_graph'].nodes(data=True)

M.analyze_reactions('basic')
M.analyze_reactions('inter')
M.analyze_reactions('up_inter')

M.total_spectrum_intensity
small_graph = M.res[7]['small_graph']

t_e_i = 0.0
for R in M.res:
    small_graph = R['small_graph']
    for G, G_D in small_graph.nodes(data=True):
        if G_D['type'] == 'G':
            print G_D['mz'].begin, G_D['mz'].end
            print G_D['intensity'], G_D['estimate']
            print
            t_e_i += G_D['estimate']
    print
    print

small_graph.edges(data=True)
