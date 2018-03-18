from __future__ import absolute_import, division, print_function
from cvxopt import matrix, spmatrix, sparse, spdiag, solvers, setseed
from six.moves import range

def diag(val, dim):
    """Create a spars matrix with val on the diagonal and zeros elsewhere."""
    return spdiag([spmatrix(val,[0],[0]) for i in range(dim)])

def normalize_rows(M):
    '''Divide rows of a matrix by their sums.'''
    for i in range(M.size[0]):
        row_hopefully = M[i,:]
        sum_abs_inv = 1.0/sum(abs(row_hopefully))
        M[i,:] = row_hopefully * sum_abs_inv
    return M
