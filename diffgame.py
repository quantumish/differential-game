import matplotlib.pyplot as plt
import numpy as np
from random import random
from numba import jit

l = 1

# @jit(nopython=True)
def euler(y_0, deriv, x_range, step=0.001):
    eps = 1e-3
    ys = [y_0]
    xs = [x_range[0]]
    # print(xs, l, f(xs[-1]), l*ys[-1])
    for i in range(int((x_range[1] - x_range[0])/step)):
        # if abs(ys[-1]) == 0: continue
        ys.append(ys[-1] + deriv(xs[-1],ys[-1])*step)
        xs.append(xs[-1] + step)
    return xs, ys

def alice_lipschitz_strat():
    return random(), 0, 0, (0, random()), (0, random())

def bob_naive_strat():
    return lambda x: x
    

for _i in range(100):
    a,b,c,p1,p2 = alice_lipschitz_strat()
    f = bob_naive_strat()
    deriv = lambda x,y: (f(x)-l*y)/(a+b*x+c*x**2)
    try:
        xs, ys = euler(p1[1], deriv, (p1[0], 1))    
        plt.plot(xs, ys)
        plt.show()
    except:
        pass
    
    # plt.scatter([params[3][0], params[4][0]], [params[3][1], params[4][1]])
    # for p in zip(xs, ys):
    #     if (abs(p[0] - params[4][0]) < 0.01) and (abs(p[1] - params[4][1]) < 0.01):
    #         print("Bob succeeds")
