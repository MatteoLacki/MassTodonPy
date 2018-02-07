import pkg_resources
import json

from MassTodonPy.Formula.Formula import Formula


def get_amino_acids(raw=True):
    """
    Retrieve the information on amino acidic bricks.

    Parameters
    ==========
    raw : boolean
        Load data from amico_acids.py?

    Returns
    =======
    A dictionary with keys (element, backbone_atom_group).
    The values are linear counters storing atom counts.
    Possible backbone_atom_group include: N, C_carbo, C_alpha.
    """
    if raw:
        from MassTodonPy.Data.amino_acids import amino_acids as amino_acids_raw
    else:
        path = pkg_resources.resource_filename('MassTodonPy', 'Data/')
        with open(path+"amino_acids.json", "rb") as f:
            amino_acids_raw = json.load(f)
    amino_acids = {
        (aa_name, brick): Formula({a: c for a, c in atom_cnt})
        for aa_name, bricks in amino_acids_raw
        for brick, atom_cnt in bricks}
    return amino_acids
