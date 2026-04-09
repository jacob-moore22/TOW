# Chapter 2: Linear Algebra

Solution of linear algebraic equations, matrix inversion, and matrix decomposition.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `gaussj` | 2.1 | Gauss-Jordan matrix inversion and linear equation solution |
| `ludcmp` | 2.3 | Linear equation solution by LU decomposition |
| `lubksb` | 2.3 | Linear equation solution by LU backsubstitution |
| `tridag` | 2.4 | Solution of tridiagonal systems |
| `mprove` | 2.5 | Iterative improvement of a linear equation solution |
| `svbksb` | 2.6 | Singular value backsubstitution |
| `svdcmp` | 2.6 | Singular value decomposition of a matrix |
| `vander` | 2.8 | Solve Vandermonde systems |
| `toeplz` | 2.8 | Solve Toeplitz systems |
| `sparse` | 2.7 | Sparse linear systems |

## Notes

- `ludcmp` + `lubksb` are used as a pair: decompose then backsubstitute.
- `svdcmp` + `svbksb` are the SVD equivalent pair.
- `mprove` refines an LU solution iteratively to improve accuracy.
- `gaussj` performs full Gauss-Jordan elimination with pivoting.
