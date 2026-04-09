# Chapter 19: Partial Differential Equations

Iterative solvers for elliptic PDEs on regular grids.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `sor` | 19.5 | Elliptic PDE solved by successive overrelaxation (SOR) method |
| `adi` | 19.6 | Elliptic PDE solved by the alternating direction implicit (ADI) method |

## Notes

- Both methods solve the same class of problem: second-order elliptic PDEs on a rectangular grid (e.g., Laplace/Poisson equations with Dirichlet boundary conditions).
- **SOR** iterates point-by-point with an overrelaxation parameter (spectral radius of the Jacobi iteration) to accelerate convergence.
- **ADI** alternates between solving implicitly in the x-direction and the y-direction each half-step, achieving faster convergence than SOR for many problems.
- The demos set up an 11x11 grid with a point source at the center and verify the solution satisfies the finite-difference equations.
