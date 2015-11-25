import numpy as np
from scipy.optimize import minimize

def func(x, sign=1.0):
	""" Objective function """
	return sign*(2*x[0]*x[1] + 2*x[0] - x[0]**2 - 2*x[1]**2)

def func_deriv(x, sign=1.0):
	""" Derivative of objective function """
	dfdx0 = sign*(-2*x[0] + 2*x[1] + 2)
	dfdx1 = sign*(2*x[0] - 4*x[1])
	return np.array([ dfdx0, dfdx1 ])

cons = (
	{	
		'type': 'eq',
		'fun' : lambda x: np.array([x[0]**3 - x[1]]),
		'jac' : lambda x: np.array([3.0*(x[0]**2.0), -1.0])
	},
	{
		'type': 'ineq',
		'fun' : lambda x: np.array([x[1] - 1]),
		'jac' : lambda x: np.array([0.0, 1.0])
	}
)

res = minimize(func, [-1.0,1.0], args=(-1.0,), jac=func_deriv, method='SLSQP', options={'disp': True})

print res

res = minimize(func, [-1.0,1.0], args=(-1.0,), jac=func_deriv, constraints=cons, method='SLSQP', options={'disp': True})
print res
