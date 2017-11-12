from collections import namedtuple

NT = namedtuple("NT", "a b")
nt = NT(a=10, b=15)



help(namedtuple)

class Blah(namedtuple, type_name="dupa" field_names=["a", "b"]):
    """Tragedy."""
