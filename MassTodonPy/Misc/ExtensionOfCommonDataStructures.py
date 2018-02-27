from collections import defaultdict

# Zondos' way of doing it.
# https://stackoverflow.com/questions/43268306/can-defaultdict-accept-callables-dependent-on-the-given-missing-key
class DefaultDict(defaultdict):
    """A defaultdict with a parametrized default."""
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError((key,))
        self[key] = value = self.default_factory(key)
        return value
