# Chapter 12: Fast Fourier Transform

FFT routines for complex, real, and multidimensional transforms, plus sine and cosine transforms.

| Program | Section | Description |
| ------- | ------- | ----------- |
| `four1`  | 12.2 | Fast Fourier transform (FFT) in one dimension |
| `twofft` | 12.3 | FFT of two real functions simultaneously |
| `realft` | 12.3 | FFT of a single real function |
| `sinft`  | 12.3 | Fast sine transform |
| `cosft`  | 12.3 | Fast cosine transform |
| `fourn`  | 12.4 | FFT in multiple dimensions |

## Notes

- `four1` is the core complex FFT -- it operates on interleaved real/imaginary data of length 2N.
- `realft` exploits symmetry to compute a real-data FFT using `four1` at half the cost.
- `twofft` computes two real FFTs simultaneously by packing them into one complex FFT.
- `sinft` and `cosft` compute discrete sine/cosine transforms via `realft`.
- `fourn` generalizes to arbitrary dimensions.
