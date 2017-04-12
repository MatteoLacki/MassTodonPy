import  numpy   as      np
from    cvxopt  import  matrix, spmatrix, sparse, spdiag, solvers

ones1 = matrix( 1.0, (3,1))
M1 = ones1 * ones1.T
ones2 = matrix( 1.0, (4,1))
M2 = ones2 * ones2.T

def make_projection(x):
    d, v = x
    ones = matrix( 1.0, (d,1))
    return v * ones * ones.T

M = matrix(spdiag(map(make_projection, [(2,5.),(3,2.),(4,10.)])))

singularValues = np.linalg.svd(np.matrix(M), compute_uv=False)
singularValues

norm(M,ord=2)

SM = spdiag(map(make_projection, [(2,5.),(3,2.),(4,10.)]))
