import pkg_resources
import json
from linearCounter.linearCounter import linearCounter as lCnt


def get_amino_acids_raw():
    """
    Retrieve the information on amino acidic bricks.

    Returns
    =======
    res: a list
    """
    path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
    with open(path+"amino_acids.json", "rb") as f:
        amino_acids = json.load(f)
    return amino_acids


def get_amino_acids():
    """
    Retrieve the information on amino acidic bricks.

    Returns
    =======
    A dictionary with keys (element, backbone_atom_group).
    The values are linear counters storing atom counts.
    Possible backbone_atom_group include: N, C_carbo, C_alpha.
    """
    amino_acids_raw = get_amino_acids_raw()
    amino_acids = {
        (aa_name, brick): lCnt({a: c for a, c in atom_cnt})
        for aa_name, bricks in amino_acids_raw
        for brick, atom_cnt in bricks}
    return amino_acids
