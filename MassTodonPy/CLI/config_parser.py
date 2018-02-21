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

# from collections import defaultdict
# import os
# from pprint import pprint
# import re
#
# from MassTodonPy.Data.Constants import eps
#
# def parse_str(name):
#     return tuple(name)
#
# def parse_float(name):
#     return name[0], float(name[1])
#
# def parse_int(name):
#     return name[0], int(name[1])
#
# def parse_modification(mod):
#     mod = mod[1:]
#     mod.reverse()
#     # Get the number of the modified amino acid.
#     AA_No = int(mod.pop().split('=')[1]) - 1 # -1: Python counts from 0
#     # Get the group to be modified.
#     group = mod.pop().split('=')[1]
#     # Get the atom diffs
#     change = {}
#     for info in mod:
#         element, diff = info.split('=')
#         change[element] = int(diff)
#     return (AA_No, group), change
#
# def parse_fragmentation_type(frags):
#     pass
#
# def parse_blocked_fragment(frags):
#     pass
#
# def parse_config_file(path):
#     '''Parse a text file in a format similar to:
#
#     * name substanceP
#     * fasta RPKPQQFFGLM
#     * charge 3
#     * modify AA=11 group=C_carbo H=1 O=-1 N=1
#     * intensity_cut_off 100
#     * mz_tol .05
#
#     The idea is actually trivial:
#         we treat each first string before space as
#         the name of the variable the value of which we want to set.
#     '''
#
#     # Vars passed to MassTodon
#     defaults = {'mz_tol': '',   # phoney
#                 'fasta': '',    # phoney
#                 'charge': 0,    # phoney
#                 'mz_digits':-10,
#                 'name':'',
#                 # 'modifications':{},
#                 # 'fragmentation_type':"cz",
#                 # 'blocked_fragments':set(['c0']),
#                 'distance_charges':5,
#                 'min_intensity':0.0,
#                 'percent_top_peaks':1.0,
#                 'deconvolution_method': 'Matteo',
#                 'joint_probability':.999,
#                 'min_prob_per_molecule':.7,
#                 '_max_buffer_len':0.5,
#                 '_L1_flow':0.001,
#                 '_L2_flow':0.001,
#                 '_L1_intensity':0.001,
#                 '_L2_intensity':0.001,
#                 '_max_times':10,
#                 '_show_progress':False,
#                 '_maxiters':1000,
#                 '_devel':False,
#                 '_sigma2':.1,
#                 '_ni2':.1}
#
#
#     types = { k: type(defaults[k]) for k in defaults }
#     complicated = ('modification', 'fragmentation_type', 'blocked_fragment')
#     for k in complicated:
#         types[k] = k
#
#     parser = {int: parse_int,
#               float: parse_float,
#               str: parse_string,
#               'modification': parse_modification,
#               'fragmentation_type': parse_fragmentation_type,
#               'blocked_fragment': parse_blocked_fragment}
#
#     with open(path, 'r') as f:
#         for line in f:
#             line = re.sub('[\n]', '', line)
#             if line:
#                 line = line.split(" ")
#                 key = line[0]
#                 if key in complicated:
#                     pass
#                 else:
#                     value = line[1]
#
#
#                 if key=='modification':
#                     AA_No, group = name
#                     args['modification']['value'][AA_No][group] = value
#                 else:
#                     args[key]['value'] = value
#
#     os.system('cls' if os.name == 'nt' else 'clear') #TODO: remove later on.
#     pprint(args)
#
#     containers = defaultdict(dict)
#     for key in args:
#         if key in ('mz_tol', ):
#             containers[key] = args[key]['value']
#         else:
#             containers[ args[key]['id'] ][key] = args[key]['value']
#
#     print()
#     pprint(containers)

    # key, value = re.sub('[\s+]|[\n]', '', line).split('=')
    # config_tmp[key] = value


    # config_tmp = {}
    # boolean_keys = ('for_plot', 'highcharts', 'raw_data', 'analyze_raw_data', 'verbose', 'csv')
    # float_keys = ('mz_prec', 'cut_off', 'opt_P', 'joint_probability_of_envelope',
    #     'min_prob_of_envelope_in_picking', 'L1_x', 'L2_x', 'L1_alpha', 'L2_alpha', 'output_deconvolution_threshold' )
    # int_keys = ('precursor_charge', 'multiprocesses_No', 'max_times_solve')
    # string_keys = ('fasta', 'frag_type', 'solver', 'method')
    #
    # with open(path, 'r') as f:
    #     for line in f:
    #         key, value = re.sub('[\s+]|[\n]', '', line).split('=')
    #         config_tmp[key] = value
    # config = {}
    # modifs = {}
    # for k in config_tmp:
    #     k_low = k.lower()
    #     if 'modification' == k_low[0:12]:
    #         aa_mod = k[12:]
    #         if not aa_mod in modifs:
    #             modifs[aa_mod] = parse_single_modification(config_tmp[k])
    #         else:
    #             print('Watch out! There are at least two modifications on the same atom:', aa_mod)
    #     elif k in boolean_keys:
    #         config[k] = parse_true_false(k, config_tmp[k])
    #     elif k in float_keys:
    #         config[k] = float(config_tmp[k])
    #     elif k in string_keys:
    #         config[k] = config_tmp[k]
    #     elif k in int_keys:
    #         config[k] = int(config_tmp[k])
    #     else:
    #         print('Key', k, 'not acceptable.')
    #         raise ValueError
    # config['modifications'] = modifs
    # return config
