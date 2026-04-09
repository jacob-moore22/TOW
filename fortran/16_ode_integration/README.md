# Chapter 16: Integration of Ordinary Differential Equations

Runge-Kutta, Bulirsch-Stoer, and adaptive step-size ODE integrators.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `rk4`    | 16.1 | Integrate one step of ODEs by fourth-order Runge-Kutta |
| `rkdumb` | 16.1 | Integrate ODEs by fourth-order Runge-Kutta (fixed step) |
| `odeint` | 16.2 | Integrate ODEs with adaptive accuracy monitoring |
| `mmid`   | 16.3 | Integrate ODEs by the modified midpoint method |
| `bsstep` | 16.4 | Integrate ODEs by the Bulirsch-Stoer method (one step) |
| `pzextr` | 16.4 | Polynomial extrapolation, used by `bsstep` |
| `rzextr` | 16.4 | Rational function extrapolation, used by `bsstep` |
| `rkqc`   | 16.2 | Runge-Kutta with quality-controlled step size |

## Notes

- `rk4` is a single fixed-step RK4 integrator. `rkdumb` calls it repeatedly with uniform steps.
- `rkqc` adds adaptive step-size control to RK4 by comparing 4th and 5th order estimates.
- `odeint` is the driver that manages step-size adaptation, calling `rkqc` or `bsstep` for each step.
- `bsstep` (Bulirsch-Stoer) uses `mmid` for modified midpoint steps and `pzextr`/`rzextr` for Richardson extrapolation. It is the recommended method for smooth problems requiring high accuracy.
