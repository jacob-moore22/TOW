# Chapter 5: Evaluation of Functions

Series summation, polynomial operations, and Chebyshev approximation.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `eulsum` | 5.1 | Sum a series by the Euler-van Wijngaarden algorithm |
| `ddpoly` | 5.3 | Evaluate a polynomial and its derivatives |
| `poldiv` | 5.3 | Divide one polynomial by another |
| `chebft` | 5.8 | Fit a Chebyshev polynomial to a function |
| `chebev` | 5.8 | Chebyshev polynomial evaluation |
| `chder`  | 5.9 | Derivative of a function already Chebyshev fitted |
| `chint`  | 5.9 | Integrate a function already Chebyshev fitted |
| `chebpc` | 5.10 | Polynomial coefficients from a Chebyshev fit |
| `pcshft` | 5.10 | Polynomial coefficients of a shifted polynomial |

## Notes

- The Chebyshev routines form a pipeline: `chebft` fits, then `chebev` evaluates, `chder` differentiates, and `chint` integrates.
- `chebpc` + `pcshft` convert from Chebyshev representation to ordinary polynomial coefficients on a shifted interval.
- `eulsum` accelerates convergence of alternating series.
