import cvxopt
import imp
import os

class cvxopt_wrapper(object):
    """ A class that should set the number of threads used by BLAS to 1.

    Notes
    -----
    This is important when running CVXOPT in any script that uses multiprocessing.
    """
    def __enter__(self):
        self.old = os.environ.get('OMP_NUM_THREADS', None)
        os.environ['OMP_NUM_THREADS'] = "1"
        imp.reload(cvxopt)

    def __exit__(self, type, value, traceback):
        if self.old:
            os.environ['OMP_NUM_THREADS'] = self.old
