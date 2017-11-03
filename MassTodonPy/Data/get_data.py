import pkg_resources
import json
import numpy as np

def get_data(dataset):
    assert dataset in ['substanceP', 'ubiquitin']
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open( path+dataset+".json","rb") as f:
        data = json.load(f)
    data["spectrum"] = tuple( np.array(d) for d in data["spectrum"] )
    return data
