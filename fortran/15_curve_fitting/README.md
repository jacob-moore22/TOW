# Chapter 15: Modeling of Data (Curve Fitting)

Linear and nonlinear least-squares fitting, robust fitting, and basis-function fitting.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `fit`    | 15.2 | Least-squares fit data to a straight line (with or without weights) |
| `lfit`   | 15.4 | General linear least-squares fit by normal equations |
| `covsrt` | 15.4 | Rearrange covariance matrix, used by `lfit` |
| `svdfit` | 15.4 | Linear least-squares fit by singular value decomposition |
| `svdvar` | 15.4 | Compute parameter variances from SVD |
| `fpoly`  | 15.4 | Polynomial basis functions for use with `lfit` or `svdfit` |
| `fleg`   | 15.4 | Legendre polynomial basis functions for use with `lfit` or `svdfit` |
| `mrqmin` | 15.5 | Nonlinear least-squares fit by Levenberg-Marquardt method |
| `mrqcof` | 15.5 | Evaluate coefficients for `mrqmin` |
| `fgauss` | 15.5 | Fit a sum of Gaussians using `mrqmin` |
| `medfit` | 15.7 | Robust straight line fit by least absolute deviation |
| `rofunc` | 15.7 | Evaluate the robust objective function, used by `medfit` |

## Notes

- `fit` is the simplest: weighted or unweighted straight-line regression.
- `lfit` and `svdfit` fit arbitrary linear basis functions (via `fpoly`, `fleg`, or user-supplied functions). `svdfit` is more numerically stable.
- `mrqmin` (Levenberg-Marquardt) handles nonlinear models -- `fgauss` demonstrates fitting a sum of Gaussians.
- `medfit` minimizes absolute deviation instead of squared deviation, making it robust to outliers.
