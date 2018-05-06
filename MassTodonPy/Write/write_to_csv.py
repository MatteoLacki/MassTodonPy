# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 Mateusz Krzysztof Łącki and Michał Startek.
#
#   This file is part of MassTodon.
#
#   MassTodon is free software: you can redistribute it and/or modify
#   it under the terms of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3.
#
#   MassTodon is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#   You should have received a copy of the GNU AFFERO GENERAL PUBLIC LICENSE
#   Version 3 along with MassTodon.  If not, see
#   <https://www.gnu.org/licenses/agpl-3.0.en.html>.


# from pandas import DataFrame as DF
# import os



# def check_and_create_folder(path):
#     """Check the path and create it if it does not exist.
#
#     Parameters
#     ----------
#     path : str
#         A path to the output directory.
#
#     Returns
#     -------
#     path : str
#         Correct path to existing directory.
#     """
#     if path[-1] != '/':
#         path += '/'
#     if not os.path.exists(path):
#         os.makedirs(path)
#     return path
#
#
#
# def i_flatten_raw_estimates(raw_estimates, threshold=0.0):
#     """Prepare a dictionary corresponding to one row in the output CSV file containing raw estimates.
#
#     Parameters
#     ----------
#     raw_estimates : list
#         A list of MassTodon raw deconvolution outputs
#     threshold : float
#         Show outputs with estimates higher than that number.
#
#     Returns
#     -------
#         out : generator
#         A generator producing dictionaries with fields: molType, formula, q, g, estimate, and status.
#     """
#     for r in raw_estimates:
#         for e in r['alphas']:
#             if e['estimate'] > threshold:
#                 e = e.copy()
#                 del e['cnt']
#                 del e['type']
#                 e['status'] = r['status']
#                 yield e
#
#
#
# def write_raw_to_csv(raw_estimates, path, threshold=0.0):
#     """Write raw deconvolution estimates to a CSV file.
#
#     Parameters
#     ----------
#     raw_estimates : list
#         A list of MassTodon raw deconvolution outputs
#     path : str
#         A path to the output directory.
#     threshold : float
#         Show outputs with estimates higher than that number.
#     """
#     try:
#         path = check_and_create_folder(path)
#         res_df = DF(i_flatten_raw_estimates(raw_estimates, threshold))
#         res_df = res_df[['molType','formula','q','g','estimate','status']]
#         res_df = res_df.sort_values('estimate', ascending=False)
#         res_df.to_csv(path_or_buf = path+'raw_estimates.csv', index = False)
#         print 'Saved estimates to file.'
#     except Exception as e:
#         print e
#         print 'Results not saved.'
#
#
#
# def iterate_algo(res, fasta):
#     """Prepare a dictionary corresponding to one row in the output CSV file containing counts and probabilities.
#
#     Parameters
#     ----------
#     res : dict
#         A dictionary of MassTodon results.
#     fasta : str
#         A string containing the amino acid sequence of the studied compound.
#
#     Returns
#     -------
#     out : generator
#         A generator producing dictionaries with fields: Algorithm, Key (description of the output), Probability, and Count.
#     """
#     for algo in ('basic_analysis','intermediate_analysis','advanced_analysis'):
#         probs, counts = res[algo]
#         for k in set(probs) | set(counts):
#             try:
#                 aa_no = int(k)
#                 aa = fasta[aa_no]
#                 L  = len(fasta)
#                 key = 'c'+str(k)+'-z'+str(L-k)+' bond broken on '+str(aa)
#             except ValueError:
#                 key = k
#
#             yield { 'Algorithm':     algo.split('_')[0],
#                     'Key':          key,
#                     'Probability':  probs.get(k, None),
#                     'Count':        counts.get(k, None)
#                    }
#
#
#
# def write_counts_n_probs_to_csv(res, fasta, path):
#     """Write counts and probabilities estimates to a CSV file.
#
#     Parameters
#     ----------
#     res : dict
#         A dictionary of MassTodon results.
#     fasta : str
#         A string containing the amino acid sequence of the studied compound.
#     path : str
#         A path to the output directory.
#     """
#     try:
#         path    = check_and_create_folder(path)
#         res_df  = DF(iterate_algo(res, fasta))[['Algorithm', 'Key', 'Count', 'Probability']]
#         res_df.to_csv(path_or_buf = path+'counts_n_probs.csv', index = False)
#         print "Saved counts'n'probs to file."
#     except Exception as e:
#         print e
#         print 'Counts not saved.'
#
#
#
# def iterate_summary(summary):
#     """Prepare a dictionary corresponding to one row in the output CSV file containing counts and probabilities.
#
#     Parameters
#     ----------
#     summary : dict
#         A summary dictionary from the standard output of the MassTodonize function, or a field (.summary) in the MassTodon class.
#
#     Returns
#     -------
#     out : generator
#         A generator producing dictionaries with fields: Statistic (e.g. L1 error), and Value (the value of the error).
#     """
#     for k in summary:
#         if not k in ('L1_error_value_error/original_total_intensity', 'L1_error_value_error/intensity_within_tolerance'):
#             yield {'Statistic': k, 'Value': summary[k]}
#
#
#
# def write_summary_to_csv(res, path):
#     """Write the summary to a CSV file.
#
#     Parameters
#     ----------
#     res : dict
#         A dictionary of MassTodon results.
#     path : str
#         A path to the output directory.
#     """
#     try:
#         path    = check_and_create_folder(path)
#         res_df  = DF(iterate_summary(res['summary']))[['Statistic', 'Value']]
#         res_df.to_csv(path_or_buf = path+'summary.csv', index = False)
#         print "Saved summary to file."
#     except Exception as e:
#         print e
#         print 'Summary not saved.'
