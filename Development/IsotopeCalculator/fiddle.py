import numpy as np
import json

try:
    import cPickle as pickle
except ImportError:
    import pickle

path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"

with open(path+"substanceP.json", "rb") as f:
    mol = json.load(f)

with open(path+"amino_acids.pickle", "rb") as f:
    bricks = pickle.load(f)

with open(path + "isotopes.pickle", "rb") as f:
    iso_masses, iso_probs = pickle.load(f)


def get_mean_and_variance(X, weights):
    X = np.array(X)
    probs = np.array(weights)
    probs = probs/sum(probs)
    average = np.dot(X, probs)
    variance = np.dot((X-average)**2, probs)
    return average, variance


means_and_variances = {el: get_mean_and_variance(iso_masses[el], iso_probs[el])
                       for el in iso_probs}


means = dict(
    (el,
     sum(pr*m for pr, m in zip(iso_probs[el],
                               iso_masses[el]))) for el in iso_masses.keys())

variances = dict(
    (el, sum( pr*m**2 for pr, m in zip(iso_probs[el], iso_masses[el])) - means[el]**2)
    for el in iso_masses.keys())


eps = 0.0000000001
for el in means_and_variances:
    mean, var = means_and_variances[el]
    assert abs(mean - means[el]) < eps
    assert abs(var - variances[el]) < eps
    # print(el)
    # print(mean, means[el])
    # print(var, variances[el])
    # print()
