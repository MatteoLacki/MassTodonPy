from scipy.optimize import linprog
import numpy as np

c = [-1, 4]
A = [[-3, 1], [1, 2]]
b = [6, 4]
x0_bounds = (None, None)
x1_bounds = (-3, None)
res = linprog(c, A_ub=A, b_ub=b, bounds=(x0_bounds, x1_bounds),
              options={"disp": True})
print(res)



res2 = linprog(c, A_eq=np.matrix(A), b_eq=np.matrix(b), bounds=(x0_bounds, x1_bounds),
              options={"disp": True})
res2
