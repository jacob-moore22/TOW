# Chapter 7: Random Numbers

Uniform, non-uniform, and cryptographic random number generators.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `ran0`   | 7.1 | Random deviate by Park and Miller minimal standard |
| `ran1`   | 7.1 | Random deviate, minimal standard plus shuffle |
| `ran2`   | 7.1 | Random deviate by L'Ecuyer long period plus shuffle |
| `ran3`   | 7.1 | Random deviate by Knuth subtractive method |
| `expdev` | 7.2 | Exponential random deviates |
| `gasdev` | 7.2 | Normally (Gaussian) distributed random deviates |
| `gamdev` | 7.3 | Gamma-law distribution random deviates |
| `poidev` | 7.3 | Poisson distributed random deviates |
| `bnldev` | 7.3 | Binomial distributed random deviates |
| `irbit1` | 7.4 | Random bit sequence |
| `irbit2` | 7.4 | Random bit sequence (alternative method) |
| `ran4`   | 7.5 | Random deviates from DES-like hashing |
| `des`    | 7.5 | Data Encryption Standard (DES) routines |
| `desks`  | 7.5 | DES key schedule |

## Notes

- `ran1` is the recommended general-purpose generator for most applications.
- `ran2` has a very long period and is suitable for large simulations.
- `gasdev` (Gaussian), `expdev` (exponential), `gamdev` (gamma), `poidev` (Poisson), and `bnldev` (binomial) all transform uniform deviates from `ran1` or `ran3`.
- `ran4` uses DES encryption (`des` + `desks`) to produce random deviates via hashing.
