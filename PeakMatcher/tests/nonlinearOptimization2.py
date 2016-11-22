import numpy as np
from scipy.optimize import minimize

def norm(X, sign=1.0):
	return sign* sum( x for x in X )

def norm_deriv(X, sign=1.0):
	return np.array([ 1., 1. ])

cons = (
	{	
		'type': 'eq',
		'fun': lambda X: np.array([ sum( x**2 for x in X) - 1 ]),
		'jac': lambda X: np.array([ 2*X[0], 2*X[1] ])
	}
)

res = minimize(fun=norm, x0=[1.0,1.0], jac=norm_deriv, constraints=cons, method='SLSQP', options={'disp': True})

print res