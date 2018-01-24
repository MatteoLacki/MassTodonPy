from collections import namedtuple
from functools import partial
from itertools import chain


def namedtuple_with_defaults(name, *fields, **fields_with_defaults):
    """Return a namedtuple with defaults set.

    We assume that 'fields' are before 'fields_with_defaults' in the
    underlying tupple, so that we can still call the factory function
    without specifying the parameters.
    """
    f = tuple(chain(fields, fields_with_defaults))
    return partial(namedtuple(name, f), **fields_with_defaults)
