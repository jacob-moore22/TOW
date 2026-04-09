# Chapter 10: Minimization and Maximization of Functions

One-dimensional minimization, multidimensional optimization, linear programming, and simulated annealing.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `mnbrak` | 10.1 | Bracket the minimum of a function |
| `golden` | 10.1 | Find minimum of a function by golden section search |
| `brent`  | 10.2 | Find minimum of a function by Brent's method |
| `dbrent` | 10.3 | Find minimum of a function using derivative information |
| `amoeba` | 10.4 | Minimize in N-dimensions by downhill simplex (Nelder-Mead) |
| `powell` | 10.5 | Minimize in N-dimensions by Powell's method |
| `linmin` | 10.5 | Minimum of a function along a ray in N-dimensions |
| `f1dim`  | 10.5 | Auxiliary function used by `linmin` |
| `frprmn` | 10.6 | Minimize in N-dimensions by conjugate gradient (Fletcher-Reeves-Polak-Ribiere) |
| `df1dim` | 10.6 | Auxiliary function (with derivatives) used by `linmin` |
| `dfpmin` | 10.7 | Minimize in N-dimensions by variable metric (BFGS) method |
| `simplx` | 10.8 | Linear programming: maximize a linear function subject to constraints |
| `simp1`  | 10.8 | Helper routine used by `simplx` |
| `simp2`  | 10.8 | Helper routine used by `simplx` |
| `simp3`  | 10.8 | Helper routine used by `simplx` |
| `anneal` | 10.9 | Traveling salesman problem by simulated annealing |
| `link`   | 10.9 | Helper routines (revcst, revers, trncst, trnspt, metrop) used by `anneal` |

## Notes

- For 1D minimization: `mnbrak` brackets, then `brent` (no derivatives) or `dbrent` (with derivatives) finds the minimum.
- `amoeba` (Nelder-Mead simplex) requires no derivatives and works well for rough landscapes.
- `frprmn` (conjugate gradient) and `dfpmin` (quasi-Newton/BFGS) require gradient information but converge faster.
- `powell` requires no derivatives and uses successive line minimizations.
- `simplx` implements the simplex method for linear programming (not related to `amoeba`'s simplex).
