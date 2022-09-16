import matplotlib.pyplot as plt
import numpy as np
from random import random
from numba import jit
from sympy import *
from zope.interface import *
from typing import Tuple, Callable, NewType, Optional

# Utility funcs
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

# globals 

l = 1

# interfaces


Point = NewType("Point", Tuple[float, float])
class AliceStrat(Interface):
    def round1() -> Tuple[float, float, float, Point, Point]:
        """Picks alpha, beta, gamma constants and the two points solutions can/must go through."""

    def round3(f: Callable[[float], float]) -> Optional[Callable[[float], float]]:
        """Picks a unique solution for the differential equation."""

class BobStrat(Interface):
    def round2(a: float, b: float, c: float, p1: Point, p2: Point) -> Optional[Tuple[Callable[[float], float], Callable[[float], float]]]:
        """Chooses a real-valued function f(x) and finds a solution to the differential equation."""

# strats

@implementer(AliceStrat)
class LipschitzAndPray():    
    def round1():
        """Hope that y(x) is not Lipschitz continuous at 0."""
        return random(), 0, 0, (0, random()), (0, random())

    def round3(f):
        """Pray it doesn't get this far - if it does, give up."""
        return None


@implementer(BobStrat)
class Analytic():
    def round2(a,b,c,p1,p2):    
        f = lambda x: 0    
        eq = Eq(Derivative(y(x), x), (x-l*y(x)/(a+b*x+c)))
        print(latex(dsolve(eq)))
        
# verifying code

def verify_diffeq(y, deriv):
    h = 0.00001
    for x in np.arange(-100, 100): # FIXME, super naive
        if (y(x+h) - y(x))/h != deriv(x,y):
            return False
    return True

alice = LipschitzAndPray()
bob = Analytic()

num_games = 10
bob_wins = [True] * num_games
for i in range(num_games):
    a,b,c,p1,p2 = alice.round1()
    f, y = bob.round2()

    deriv = lambda x,y: (f(x)-l*y(x))/(a+b*x+c**2)
    if verify_diffeq(y, deriv):    
        y2 = alice.round3(f)        
        if verify_diffeq(y2, deriv):
            bob_wins[i] = False
    else:
        bob_wins[i] = False

print(bob_wins)
