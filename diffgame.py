import matplotlib.pyplot as plt
import numpy as np
from random import *
import math
seed(10)
from numba import jit
from sympy import *
from zope.interface import *
from typing import Tuple, Callable, NewType, Optional
from tqdm import tqdm
from enum import Enum
from matplotlib.animation import FuncAnimation
import sys

# Utility funcs
def euler(y_0, deriv, x_range, step=0.001):
    eps = 1e-3
    ys = [y_0]
    xs = [x_range[0]]
    # print(xs, l, f(xs[-1]), l*ys[-1])
    for i in range(int((x_range[1] - x_range[0]) / step)):
        # if abs(ys[-1]) == 0: continue
        ys.append(ys[-1] + deriv(xs[-1], ys[-1]) * step)
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
    def round2(
        a: float, b: float, c: float, p1: Point, p2: Point
    ) -> Optional[Tuple[Callable[[float], float], Callable[[float], float]]]:
        """Chooses a real-valued function f(x) and finds a solution to the differential equation."""


# strats


@implementer(AliceStrat)
class LipschitzAndPray:
    def round1():
        """Hope that y(x) is not Lipschitz continuous at 0."""
        return 0, random(), 0, (0, random()), (0, random())

    def round3(f):
        """Pray it doesn't get this far - if it does, give up."""
        return None


@implementer(AliceStrat)
class NaiveRoots:
    def round1():
        root = random()
        b, a = uniform(0, root), uniform(root, 1)
        return root**2, -2*root, 1, (b, random()), (a, random())

    def round3(f):
        return None

@implementer(AliceStrat)
class Roots:
    def round1():
        root1, root2 = random(), random()
        while root1 == root2: # reinit if unlucky
            root1, root2 = random(), random()            
        return root1*root2, -(root1+root2), 1, (root1, random()), (root2, random())

    def round3(f):
        return None


@implementer(AliceStrat)
class Random:
    def round1():
        return random(), random(), random(), (random(), random()), (random(), random())

    def round3(f):
        return None

@implementer(AliceStrat)
class Rude:
    """Be mean."""

    def round1():
        big = sys.float_info.max
        return big, big, big, (big, big), (big, big)

    def round3(f):
        return sys.float_info.min


@implementer(BobStrat)
class PrecomputedAnalyticZero:
    """Assumes a,b,c are nonzero and picks f(x) = 0 to make things simple.
    Plugs into a precomputed analytical solution to get y(x).
    """
    def round2(a, b, c, p1, p2):
        try:
            scary_term = lambda x: math.exp(
                -(2 * l * math.atan((b + 2 * c * x) / math.sqrt(4 * a * c - b**2)))
                / math.sqrt(4 * a * c - b**2)
            )
            k = p1[1] / scary_term(p1[0])
            return lambda x: 0, lambda x: k * scary_term(x)
        except:
            return None


@implementer(BobStrat)
class Analytic:
    """Uses SymPy to solve the differential equation."""
    def __init__(self, f):
        self.f = f

    def round2(self, a, b, c, p1, p2):
        y = Function("y")
        sa, sb, sc, sl, sx = symbols("a b c l x")
        eq = Eq(Derivative(y(sx), sx), (self.f-sl * y(sx) / (sa + sb * sx + sc * sx**2)))
        eq = eq.subs({sa: a, sb: b, sc: c, sl: l})
        sol = dsolve(eq)
        k = p1[1]/sol.subs({Symbol("C1"): 1, sx: p1[0]}).rhs
        sol = sol.subs({Symbol("C1"): k})
        return lambda x: 0, lambda x: sol.subs(sx, x).rhs


@implementer(BobStrat)
class Euler:
    def __init__(self, f):
        self.f = f

    """Approximates a solution numerically using Euler's method."""
    def round2(self,a,b,c,p1,p2):
        try:
            deriv = lambda x, y: (self.f(x) - l * y) / (a + b * x + c * x**2)
            xs, ys = euler(p1[1], deriv, (p1[0], 1))
        except:
            return None

        # print(ys)
        def near_analytic(x):
            for i in range(len(xs))[:-1]:
                if x >= xs[i] and x < xs[i+1]:
                    return ys[i]
            return np.nan

        return self.f, near_analytic

# verifying code

def verify_diffeq(y, deriv, p1, p2):
    """Verifies a potential solution.
    - Checks that it satisfies the differential equation
    - Checks that it passes through at least one of the points
    """
    # verify it satisfies the diffeq
    h = 0.001
    for x in np.arange(p1[0], p2[0], 0.01): # FIXME, super naive
        if abs((y(x + h) - y(x)) / h - deriv(x, y)) > 0.001:
            # print((y(x+h) - y(x))/h, deriv(x,y))
            return None

    epsilon = 0.001
    return (abs(y(p1[0])-p1[1]) < epsilon, abs(y(p2[0]) - p2[1]) < epsilon)

class Result(Enum):
    WIN = (1,)
    GIVEUP = (2,)
    NO_POINT = (3,)
    INVALID_EQ = (4,)
    ALICE_WIN = (5,)

    def __str__(self):
        return self.name

alice = Random
bob = Euler(lambda x: 0)
    

fig = plt.figure()
ax = plt.axes()
ax.set_title("Random (Alice) vs. Euler0 (Bob)")
ax.set_xlim(-.1, 1.1)
ax.set_ylim(-.1, 1.1)
line, = ax.plot([], [], lw=2)
scatter = ax.scatter([], [])
res_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
stats_text = ax.text(0.05, 0.05, '', transform=ax.transAxes)
results = []

def play_game(i):
    a, b, c, p1, p2 = alice.round1()
    if p1[0] <= 0 or p1[1] <= 0 or p2[0] <= 0 or p2[1] <= 0:
        return # TODO make less bad
    scatter.set_offsets([p1, p2])

    bob_choice = bob.round2(a, b, c, p1, p2)
    if bob_choice != None:        
        f, y = bob_choice

        result = Result.WIN

        xs = np.arange(0, 1, 0.01)
        line.set_data(xs, list(map(y, xs)))        
        deriv = lambda x, y: (f(x) - l * y(x)) / (a + b * x + c * x**2)    
        pts = verify_diffeq(y, deriv, p1, p2)    
        if pts is not None:
            # TODO implement forcing Alice to go through the same points as Bob
            if True in pts:
                y2 = alice.round3(f)
                if y2 is not None and verify_diffeq(y2, deriv, p1, p2):
                    result = Result.ALICE_WIN
            else:
                result = Result.NO_POINT
        else:
            result = Result.INVALID_EQ
    else:
        result = Result.GIVEUP
        
    results.append(result)
    res_text.set_text(str(result))
    res_text.set_c("g" if str(result) == "WIN" else "r")

    wins = round(results.count(Result.WIN)/len(results) * 100, 2)
    no_points = round(results.count(Result.NO_POINT)/len(results) * 100, 2)
    giveups = round(results.count(Result.GIVEUP)/len(results) * 100, 2)
    invalid_eqs = round(results.count(Result.INVALID_EQ)/len(results) * 100, 2)
    alice_wins = round(results.count(Result.ALICE_WIN)/len(results) * 100, 2)
    
    stats_text.set_text(f"WIN {wins}%\nNO_POINT {no_points}%\nGIVEUP {giveups}%\nINVALID_EQ {invalid_eqs}%\nALICE_WIN {alice_wins}%")
    
    return [line, scatter, res_text, stats_text]

anim = FuncAnimation(fig, play_game, frames=200, interval=20, blit=True)
anim.save("./demo.gif")
