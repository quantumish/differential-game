# The Game

1. Alice chooses three scalars, $\alpha, \beta, \text{ and } \gamma$, as well as two points, $(x_1, y_1)$ and $(x_2, y_2)$, in the first quadrant.
2. Bob chooses a function $f(x)$ that satisfies $\lim_{x \rightarrow 0} f(x) = b < \infty$ such that there exists a solution to $(\alpha + \beta x + \gamma x^2)y' + \lambda y = f(x)$ that passes through one or both of ${(x_1, y_1), (x_2, y_2)}$. Let the set of points that Bob's solution passes through be $S$.
3. Alice finds a distinct solution from Bob's that passes through the points in $S$.

The last player to make a move wins. That is, Bob can win by Alice failing at step 3, while Alice can win either by Bob failing at step 2 or by herself succeeding at step 3.

We want to determine who wins with what strategies for which values of $\lambda$ (set before the game as a parameter) and design a computational environment to let us experience the fun of differential games ourselves.

# Special Cases
## 1. $\beta = \gamma = 0, \alpha \neq 0$
### 1.1 Bob Wins
*Proposition A:* Under special case (1), Bob is guaranteed to win.

**Proof A:**
We have, after rearranging, $$y' = \frac{f(x) - \lambda y}{\alpha},$$ which follows the standard form of $y' = F(x, y)$. If Bob chooses $f(x) = x$, which satisfies the condition by $\lim_{x\rightarrow 0} x = 0 < \infty$, we have $F(x, y) = \frac{x-\lambda y}{\alpha}$. Though we could use Edwards' & Penney's Version of the Existence/Uniqueness Theorem, which is computationally easier, this is a simple enough case that we can prove the slightly more general Lipschitz Continuity and apply Taylor's Version of the Existence/Uniqueness Theorem. See 1.2 for that proof. The result allows us to apply Taylor's Version of the Existence/Uniquness Theorem. For brevity, we'll refer to this as "EUT (Taylor's Version)", as opposed to "EUT (Edwards' Version)". Therefore, there is one solution through either of the points that satisfies the original differential equation. Bob chooses this, and because the solution is unique, Alice will be unable to reply. $\mathcal{QED}\quad\blacksquare$

### 1.2 EOT (Taylor's Version)
The conditions for EOT (Taylor's Version) to apply to $y' = F(x, y)$ are as follows.

1. $F(x, y)$ is bounded and continous on $I \times \Omega$ where $I$ is an open interval containing $x_0$ and $\Omega$ is an open interval containing $y_0$, with $(x_0, y_0)$ being an initial condition.
2. $F(x, y)$ is Lipschitz Continuous in y, that is: $$\left|\frac{F(x, y) - F(x, y_0)}{y - y_0}\right| \leq L \quad \forall x \in I; y, y_0 \in \Omega; L \in (0, \infty).$$

We posited in 1.1 that Bob could satisfy EOT (Taylor's Version) by choosing $f(x) = x$.First, let's find an expression for $F(x, y)$.

$$F(x, y) = y' = \frac{f(x) - \lambda y}{\alpha} = \frac{x - \lambda y}{\alpha}$$

Obviously, this function is bounded and continuous on any finite interval, so the first condition is satisfied.

To tackle the second condition, let's evaluate the right side of the inequality above and show that $F(x, y)$ is Lipschitz Continuous for some finite $L$.

$$\left|\frac{\frac{x - \lambda y}{\alpha} - \frac{x - \lambda y_0}{\alpha}}{y - y_0}\right| = \left|\frac{\lambda}{\alpha}\frac{y - y_0}{y - y_0}\right| = \left|\frac\lambda\alpha\right|$$

Thus, we can simply choose $L = \left|\frac\lambda\alpha\right|$ as an upper bound for the growth of $F(x, y)$. With both conditions met, EOT (Taylor's Version) is satisfied. $\mathcal{QED} \quad \blacksquare$

## 2. $\alpha = \gamma = 0, \beta \neq 0$
### 2.1 Bob Wins
*Proposition B:* Under special case (2), Bob is guaranteed to win.

**Proof B:** We'll go through this proof in an expedited fashion; if unclear, supplement with information from 1.1, as the proof follows the same general structure. Using the scalars' values from the special case, we derive $F(x, y) = \frac{f(x) - \lambda y}{\beta x}$. Here, let's use EUT (Edwards' Version). We can take $D_y(F(x, y)) = -\frac{\lambda}{\beta x}$, which is defined for all $x, y \in I \times \Omega$. Thus, a unique solution exists to the differential equation. Bob chooses this, and because the solution is unique, Alice will be unable to reply. $\mathcal{QED} \quad \blacksquare$

## 3. TBD

## 4. TBD