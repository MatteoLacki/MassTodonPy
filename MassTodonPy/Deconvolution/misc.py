from    networkx.linalg.attrmatrix  import attr_matrix
import  numpy                       as     np
from    scipy.optimize              import nnls

def get_matrix_representation(G, total_intensities, min_mz, max_mz):
    """Get the representation of the problem in terms of Y and X matrices."""
    mol_columns = np.array([N < 0  for N in G])
    peak_rows   = np.array([N >= 0 for N in G])
    X, ordering = attr_matrix(G, edge_attr='prob')
    X = X[:,mol_columns][peak_rows,:]
    ordering = np.array(ordering)
    peaks    = ordering[ordering >= 0]
    Y = np.concatenate((total_intensities[peaks],
                        np.zeros(X.shape[1])))
    x = 1.0 - np.array(X.sum(axis=0)).flatten()
    X = np.concatenate((X,
                        np.diag(x)))
    return Y, X, min_mz[peaks], max_mz[peaks]