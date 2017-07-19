import cPickle as pickle
from pandas import DataFrame as DF
from MassTodonPy.Outputing.write_to_csv import write_raw_to_csv, i_flatten_raw_estimates

datapath = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/ubi_only_real/'
with open(datapath+'ubiquitins.masstodon', 'r') as h:
    results = pickle.load(h)

results[0].keys()


def iter_results(results):
    for res in results:
        settings = res['experimental_setting']
        for estimate in i_flatten_raw_estimates(res['masstodon_output']['raw_estimates']):
            settings_tmp = settings.copy()
            settings_tmp.update(estimate)
            yield settings_tmp

DF(iter_results(results)).to_csv( path_or_buf = datapath+'raw_estimates.csv', index = False )
