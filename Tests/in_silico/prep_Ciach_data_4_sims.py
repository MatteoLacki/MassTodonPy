# from deconv_misc import change_key, sigmas2probs, probs2sigmas
#
# sigmas = [probs2sigmas[a] for a in (0.01168997000000005, 0.14815520000000004, 0.49865629)]
# # fp_main= sys.argv[1]
# fp_main='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico'
#
# fp_in  = fp_main+'/results_Ciach/'
# fp_out = fp_main+'/results_Matteo/'
# molsNo = 100000
#
#
# i = 0
# for molsNo in (1000, 10000, 100000):
#     with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
#         res = pickle.load(f)
#     for r in res:
#

# def helper(helper_args):
#     helper_args
#
# def run(molsNo, fp_in, fp_out):
#
#     P = Pool(multiprocesses_No)
#     results = P.map(    helper,
#                         pool_args    )
#
#     P.close()
#     P.join()
#
#     with open(fp_in+'results_molsNo-'+str(molsNo), "rb") as f:
#         res = pickle.load(f)
#     RES = []
