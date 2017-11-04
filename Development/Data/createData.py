from MassTodonPy import MassTodonize, get_data

substanceP = get_dataset('substanceP')
ubiquitin  = get_dataset('ubiquitin')


def serialize(mol):
    mol = mol.copy()
    masses, intensities = mol["spectrum"]
    mol["spectrum"] = masses.tolist(), intensities.tolist()
    return mol

data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"
with open(data_path + "substanceP.json", "wb") as f:
    json.dump(serialize(substanceP), f)

with open(data_path + "ubiquitin.json", "wb") as f:
    json.dump(serialize(ubiquitin), f)
