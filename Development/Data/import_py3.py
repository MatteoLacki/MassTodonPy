import pickle

data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"
with open(data_path + "substanceP.example", "rb") as f:
    mol = pickle.load(f)
