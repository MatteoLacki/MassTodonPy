from  deconvolution_tests_misc import simulate, add_quants, getMols
import sys
import cPickle as pickle

sigma = float(sys.argv[1])
path_to_results = sys.argv[2]
min_x0_exp = int(sys.argv[3])
max_x0_exp = int(sys.argv[4])
maxMolsNo  = int(sys.argv[5])

mzPrec = 0.05
cutOff = 100

# sigma = 0.03
# path_to_results ='/Users/matteo/Documents/MassTodon/in_silico_results/'
# min_x0_exp, max_x0_exp = 3, 3
# maxMolsNo = 1

for molsNo in xrange(1, maxMolsNo+1):
    partial_results = []
    for x0_exp in xrange(min_x0_exp, max_x0_exp+1):
        x0 = 10 ** x0_exp
        mols = [('p','C378H629N105O118S1',76,1,i) for i in xrange(molsNo)]
        mols = add_quants(x0, mols)
        res  = simulate(mols, sigma, Q=molsNo+1)
        partial_results.append( (mols, x0, sigma, res) )

    filepath = path_to_results+'sigma_'+str(sigma)[2:]+'_molsNo_'+str(molsNo)+'.matteo'
    print filepath
    with open(filepath, 'w') as f:
        pickle.dump( partial_results, f)

# with open(path_to_results+'sigma_999999426697_molsNo_1.matteo', 'rb') as f:
#     x = pickle.load(f)


#
# x0 = 10 ** x0_exp
# mols = [('p','C378H629N105O118S1',76,1,i) for i in xrange(molsNo)]
# mols = add_quants(x0, mols)
# res  = simulate(mols, sigma, Q=molsNo+1)
# partial_results.append( (mols, x0, sigma, res) )
