import matplotlib.pyplot as plt
import numpy as np
import random
import math

random.seed(10)
from random import random
from numba import jit
from sympy import *
from zope.interface import *
from typing import Tuple, Callable, NewType, Optional
from tqdm import tqdm
from enum import Enum
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


# import traceback
@implementer(BobStrat)
class PrecomputedAnalytic:
    def round2(a, b, c, p1, p2):
        try:
            scary_term = lambda x: math.exp(
                -(2 * l * math.atan((b + 2 * c * x) / math.sqrt(4 * a * c - b**2)))
                / math.sqrt(4 * a * c - b**2)
            )
            k = p1[1] / scary_term(p1[0])
            return lambda x: 0, lambda x: k * scary_term(x)
        except:
            # traceback.print_exc()
            # print(a,b,c,l,p1[0])
            return None


@implementer(BobStrat)
class Analytic:
    def round2(a, b, c, p1, p2):
        f = lambda x: 0
        y = Function("y")
        sa, sb, sc, sl, sx = symbols("a b c l x")
        eq = Eq(Derivative(y(sx), sx), (-sl * y(sx) / (sa + sb * sx + sc * sx**2)))
        eq = eq.subs({sa: a, sb: b, sc: c, sl: l})
        sol = dsolve(eq)
        k = p1[1]/sol.subs({Symbol("C1"): 1, sx: p1[0]}).rhs
        sol = sol.subs({Symbol("C1"): k})         
        # print(sol.rhs)
        return f, lambda x: sol.subs(sx, x).rhs


# verifying code


def verify_diffeq(y, deriv, p1, p2):
    # verify it satisfies the diffeq
    h = 0.00001
    for x in np.arange(0, 1, 0.01):  # FIXME, super naive
        if abs((y(x + h) - y(x)) / h - deriv(x, y)) > 0.01:
            # print(x, y(x), (y(x+h) - y(x))/h, deriv(x,y))
            return None

    return (y(p1[0]) == p1[1], y(p2[0]) == p2[1])


class Result(Enum):
    WIN = (1,)
    GIVEUP = (2,)
    NO_POINT = (3,)
    INVALID_EQ = (4,)
    ALICE_WIN = (5,)


alice = Random
bob = Analytic

num_games = 100
bob_wins = [Result.WIN] * num_games
for i in tqdm(range(num_games)):
    a, b, c, p1, p2 = alice.round1()
    bob_choice = bob.round2(a, b, c, p1, p2)
    if bob_choice == None:
        bob_wins[i] = Result.GIVEUP
        continue
    f, y = bob_choice
    # xs = np.arange(-100, 100)
    # plt.plot(xs, list(map(y, xs)))
    # plt.show()

    deriv = lambda x, y: (f(x) - l * y(x)) / (a + b * x + c * x**2)
    pts = verify_diffeq(y, deriv, p1, p2)
    if pts == None:
        bob_wins[i] = Result.INVALID_EQ
    elif True in pts:
        y2 = alice.round3(f)
        if y2 is not None and verify_diffeq(y2, deriv, p1, p2):
            bob_wins[i] = Result.ALICE_WIN
    else:
        bob_wins[i] = Result.NO_POINT


print(f"Bob wins {bob_wins.count(Result.WIN)/num_games * 100}% of the time!")
print(f"   - lost to points {bob_wins.count(Result.NO_POINT)/num_games * 100}% of the time!")
print(f"   - lost to giveup {bob_wins.count(Result.GIVEUP)/num_games * 100}% of the time!")
print(f"   - lost to invalid {bob_wins.count(Result.INVALID_EQ)/num_games * 100}% of the time!")
print(f"   - lost to alice {bob_wins.count(Result.ALICE_WIN)/num_games * 100}% of the time!")
