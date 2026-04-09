# Chapter 6: Special Functions

Gamma, beta, error, Bessel, and elliptic functions, plus factorials and binomial coefficients.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `gammln` | 6.1 | Logarithm of the gamma function |
| `factrl` | 6.1 | Factorial function |
| `bico`   | 6.1 | Binomial coefficients |
| `factln` | 6.1 | Logarithm of the factorial function |
| `beta`   | 6.1 | Beta function |
| `gammp`  | 6.2 | Incomplete gamma function P(a,x) |
| `gammq`  | 6.2 | Complement of incomplete gamma function Q(a,x) = 1 - P(a,x) |
| `gser`   | 6.2 | Series expansion used by `gammp` and `gammq` |
| `gcf`    | 6.2 | Continued fraction used by `gammp` and `gammq` |
| `erf`    | 6.2 | Error function erf(x) |
| `erfc`   | 6.2 | Complementary error function erfc(x) = 1 - erf(x) |
| `erfcc`  | 6.2 | Complementary error function, concise routine |
| `betai`  | 6.4 | Incomplete beta function |
| `betacf` | 6.4 | Continued fraction used by `betai` |
| `bessj0` | 6.5 | Bessel function J0 |
| `bessy0` | 6.5 | Bessel function Y0 |
| `bessj1` | 6.5 | Bessel function J1 |
| `bessy1` | 6.5 | Bessel function Y1 |
| `bessy`  | 6.5 | Bessel function Y of general integer order |
| `bessj`  | 6.5 | Bessel function J of general integer order |
| `bessi0` | 6.6 | Modified Bessel function I0 |
| `bessk0` | 6.6 | Modified Bessel function K0 |
| `bessi1` | 6.6 | Modified Bessel function I1 |
| `bessk1` | 6.6 | Modified Bessel function K1 |
| `bessk`  | 6.6 | Modified Bessel function K of integer order |
| `bessi`  | 6.6 | Modified Bessel function I of integer order |
| `plgndr` | 6.8 | Associated Legendre polynomials (spherical harmonics) |
| `sncndn` | 6.11 | Jacobian elliptic functions sn, cn, dn |
| `cel`    | 6.12 | Complete elliptic integral |
| `el2`    | 6.12 | Elliptic integral of the second kind |

## Notes

- `gammln` is a foundational routine called by many others (`factrl`, `factln`, `bico`, `gammp`, `gammq`, `gser`, `gcf`, `betai`, `bnldev`, `poidev`).
- `gammp`/`gammq` use `gser` for small x and `gcf` for large x.
- `erf` and `erfc` are built on top of `gammp`/`gammq`.
- Bessel functions J and Y (ordinary) and I and K (modified) are provided for orders 0, 1, and general integer order.
