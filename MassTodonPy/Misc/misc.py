import cvxopt
import imp
import os

class cvxopt_wrapper(object):
    def __enter__(self):
        self.old = os.environ.get('OMP_NUM_THREADS', None)
        os.environ['OMP_NUM_THREADS'] = "1"
        imp.reload(cvxopt)

    def __exit__(self):
        os.environ['OMP_NUM_THREADS'] = self.old
