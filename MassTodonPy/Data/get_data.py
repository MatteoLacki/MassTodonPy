import pkg_resources
import json
import numpy as np
try:
    import cPickle as pickle
except ImportError:
    import pickle


def get_dataset(dataset):
    assert dataset in ['substanceP', 'ubiquitin']
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path+dataset+".json", "rb") as f:
        data = json.load(f)
    data["spectrum"] = tuple(np.array(d) for d in data["spectrum"])
    return data


def get_bricks():
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path+"amino_acids.txt", "rb") as f:
        bricks = pickle.load(f)
    return bricks
