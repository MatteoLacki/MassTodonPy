from six.moves import range, zip

from MassTodonPy.Data.Constants import infinity


def get_buffer(r_prev, l, r, l_next, max_length=.5):
    """Create one buffer."""
    tol = min(max_length, (l-r_prev)/2, (l_next-r)/2)
    return l - tol, r + tol

def buffer_args(L, R, max_length=.5):
    """Get arguments to create different buffers."""
    yield -infinity, L[0], R[0], L[1], max_length
    for i in range(1,len(L)-1):
        yield R[i-1], L[i], R[i] , L[i+1], max_length
    yield R[-2], L[-1], R[-1], infinity, max_length

def get_buffers(L, R, max_length=.5):
    """Get left and right ends of buffers around intervals.

    Parameters
    ==========
    L : iterable
        Left ends of intervals.
    R : iterable
        Right ends of intervals.
    max_length : float
        Maximal length of a buffer.
    Returns
    =======
    out : iterable
        Two lists of left and right ends of buffers.
    """
    return zip(*(get_buffer(*args) for args in buffer_args(L, R, max_length)))
