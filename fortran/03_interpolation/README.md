# Chapter 3: Interpolation and Extrapolation

Polynomial, rational, and spline interpolation in one and two dimensions.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `polint` | 3.1 | Polynomial interpolation (Neville's algorithm) |
| `ratint` | 3.2 | Rational function interpolation |
| `spline` | 3.3 | Construct a cubic spline |
| `splint` | 3.3 | Cubic spline interpolation (evaluate the spline) |
| `locate` | 3.4 | Search an ordered table by bisection |
| `hunt`   | 3.4 | Search a table when calls are correlated |
| `polcoe` | 3.5 | Polynomial coefficients from a table of values |
| `polcof` | 3.5 | Polynomial coefficients from a table of values (alternative) |
| `polin2` | 3.6 | Two-dimensional polynomial interpolation |
| `bcucof` | 3.6 | Construct two-dimensional bicubic coefficients |
| `bcuint` | 3.6 | Two-dimensional bicubic interpolation |
| `splie2` | 3.6 | Construct two-dimensional spline |
| `splin2` | 3.6 | Two-dimensional spline interpolation |

## Notes

- `spline` + `splint` are a pair: first construct the spline, then evaluate it.
- `splie2` + `splin2` are their two-dimensional equivalents.
- `locate` and `hunt` are table-search utilities used by many other routines.
- `bcucof` + `bcuint` handle bicubic interpolation on a grid.
