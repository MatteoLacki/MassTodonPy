import itertools


def alternate(*args):
    """Alternate between the elements from different iterators."""
    for x in itertools.chain(*zip(*args)):
        yield x
