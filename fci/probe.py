import numpy as np
from skmonaco import mcquad
# f = lambda x_y,alpha,beta: np.exp(-alpha*x_y[0]**2)*np.exp(-beta*x_y[1]**2)

def f(x_y, alpha, beta):
    return np.exp(-alpha*x_y[0]**2)*np.exp(-beta*x_y[1]**2)


alpha = 1.0
beta = 2.0
ans, err = mcquad(f,xl=[0.,0.],xu=[1.,1.],npoints=100000,args=(alpha,beta))
print ans
print err
