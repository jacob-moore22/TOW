# Chapter 9: Root Finding and Nonlinear Sets of Equations

Bracketing, bisection, Brent's method, Newton-Raphson, and polynomial root finding.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `scrsho` | 9.0 | Graph a function to visually search for roots |
| `zbrac`  | 9.1 | Outward search for brackets on roots |
| `zbrak`  | 9.1 | Inward search for brackets on roots |
| `rtbis`  | 9.1 | Find root of a function by bisection |
| `rtflsp` | 9.2 | Find root of a function by false-position (regula falsi) |
| `rtsec`  | 9.2 | Find root of a function by secant method |
| `zbrent` | 9.3 | Find root of a function by Brent's method |
| `rtnewt` | 9.4 | Find root of a function by Newton-Raphson |
| `rtsafe` | 9.4 | Find root by Newton-Raphson with bisection safeguard |
| `laguer` | 9.5 | Find a root of a polynomial by Laguerre's method |
| `zroots` | 9.5 | All roots of a polynomial by Laguerre's method with deflation |
| `qroot`  | 9.5 | Complex or double root of a polynomial by Bairstow's method |
| `mnewt`  | 9.6 | Newton's method for systems of nonlinear equations |

## Notes

- `zbrent` (Brent's method) is generally the recommended scalar root finder -- it combines bisection reliability with superlinear convergence.
- `rtsafe` is preferred when derivatives are available -- it combines Newton-Raphson speed with bisection safety.
- `zbrac`/`zbrak` are used to find bracketing intervals before calling a root finder.
- `zroots` finds all roots of a polynomial; `laguer` finds one root at a time.
- `mnewt` generalizes Newton's method to vector-valued functions (systems of equations).
