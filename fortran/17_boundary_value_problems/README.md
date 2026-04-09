# Chapter 17: Two-Point Boundary Value Problems

Shooting methods and relaxation methods for ODEs with boundary conditions at two endpoints.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `shoot`  | 17.1 | Solve a two-point boundary value problem by shooting |
| `shootf` | 17.2 | Solve by shooting to a fitting point |
| `solvde` | 17.3 | Solve a two-point BVP by relaxation on a mesh |
| `bksub`  | 17.3 | Backsubstitution, used by `solvde` |
| `pinvs`  | 17.3 | Diagonalize a sub-block, used by `solvde` |
| `red`    | 17.3 | Reduce columns of a matrix, used by `solvde` |
| `sfroid` | 17.4 | Spheroidal functions computed by the method of `solvde` (standalone program) |
| `difeq`  | 17.4 | Spheroidal matrix coefficients, used by `sfroid` |

## Notes

- **Shooting** (`shoot`, `shootf`): Convert the BVP into an initial value problem and iterate. `shoot` fires from one boundary; `shootf` fires from both and matches at an interior point.
- **Relaxation** (`solvde`): Discretize the ODE on a mesh and solve the resulting nonlinear system iteratively. `bksub`, `pinvs`, and `red` are its internal linear algebra helpers.
- `sfroid` is a standalone `PROGRAM` that demonstrates `solvde` by computing spheroidal wave functions. `difeq` supplies the ODE coefficients for this problem.
