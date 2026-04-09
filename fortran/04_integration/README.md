# Chapter 4: Integration of Functions

Numerical quadrature methods from basic trapezoidal to adaptive Romberg and Gaussian schemes.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `trapzd` | 4.2 | Trapezoidal rule (single refinement step) |
| `qtrap`  | 4.2 | Integrate using the trapezoidal rule to convergence |
| `qsimp`  | 4.2 | Integrate using Simpson's rule |
| `qromb`  | 4.3 | Integrate using Romberg adaptive method |
| `midpnt` | 4.4 | Extended midpoint rule |
| `qromo`  | 4.4 | Integrate using open Romberg adaptive method |
| `midinf` | 4.4 | Integrate a function on a semi-infinite interval |
| `qgaus`  | 4.5 | Integrate a function by Gaussian quadratures |
| `gauleg` | 4.5 | Compute Gauss-Legendre weights and abscissas |
| `quad3d` | 4.6 | Integrate a function over a three-dimensional space |

## Notes

- `trapzd` is the building block: it refines a trapezoidal estimate by one level. `qtrap`, `qsimp`, and `qromb` all call `trapzd` iteratively.
- `qromb` uses Romberg extrapolation (via `polint`) for faster convergence.
- `qromo` uses open formulas (via `midpnt` or `midinf`) for improper integrals.
- `gauleg` precomputes Gauss-Legendre nodes/weights for use with `qgaus`.
