from MassTodonPy import get_dataset
import json

substanceP = get_dataset('substanceP')


# substanceP['modifications'] = {11: {"C_carbo": substanceP['modifications'][(11,"C_carbo")]}}
# substanceP['modifications'] = {int(k): v for k, v in substanceP['modifications'].items()}
ubiquitin = get_dataset('ubiquitin')


def serialize(mol):
    mol = mol.copy()
    masses, intensities = mol["spectrum"]
    mol["spectrum"] = masses.tolist(), intensities.tolist()
    return mol

data_path = "/Users/matteo/Documents/MassTodon/MassTodonPy/MassTodonPy/Data/"
with open(data_path + "substanceP.json", "w") as f:
    json.dump(serialize(substanceP), f)

with open(data_path + "ubiquitin.json", "wb") as f:
    json.dump(serialize(ubiquitin), f)
