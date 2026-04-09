# Chapter 13: Fourier and Spectral Applications

Convolution, correlation, spectral estimation, linear prediction, and data smoothing.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `convlv` | 13.1 | Convolution or deconvolution of data using FFT |
| `correl` | 13.2 | Correlation or autocorrelation of data using FFT |
| `spctrm` | 13.4 | Power spectrum estimation using FFT (Welch's method) |
| `memcof` | 13.6 | Evaluate maximum entropy method (MEM) coefficients |
| `fixrts` | 13.6 | Reflect polynomial roots into the unit circle (stabilize MEM filter) |
| `predic` | 13.6 | Linear prediction using MEM coefficients |
| `evlmem` | 13.7 | Power spectral estimation from MEM coefficients |
| `smooft` | 13.1 | Smooth data in the frequency domain using FFT |

## Notes

- `convlv` and `correl` both use FFT to achieve O(N log N) performance.
- `memcof` + `evlmem` provide maximum entropy spectral estimation, which gives sharper peaks than FFT-based methods for short data records.
- `memcof` + `fixrts` + `predic` form a linear prediction pipeline: estimate AR coefficients, stabilize the filter, then extrapolate.
- `spctrm` uses segment averaging (Welch's method) for robust spectral estimates.
