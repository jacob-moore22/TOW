# CONVLV — FFT-Based Convolution (MATAR Parallel)

Parallel C++ conversion of the Numerical Recipes `CONVLV` example using
[MATAR](https://github.com/lanl/MATAR) with a Kokkos backend.  Demonstrates
cross-chapter dependency handling (four subroutines across two chapters) and
five distinct parallelization patterns.

## Dependency Graph

```
convlv.dem (driver)
  └─ convlv.f  (ch13 — spectral analysis)
       ├─ twofft.f (ch12 — Fourier transform)
       │    └─ four1.f (ch12)
       └─ realft.f  (ch12)
            └─ four1.f (ch12)
```

Each subroutine is converted into its own header (`four1.hpp`, `realft.hpp`,
`twofft.hpp`, `convlv.hpp`) so the cross-chapter structure is preserved.

## Algorithm

`CONVLV` performs convolution in the frequency domain:

1. Pad the response function with wrap-around symmetry and zero-fill.
2. Simultaneously FFT signal and response via `TWOFFT` → `FOUR1`.
3. Pointwise complex multiply in frequency domain.
4. Inverse real-FFT via `REALFT` → `FOUR1` to recover the convolution.

The driver validates by comparing the FFT result against direct circular
convolution.

## Parallelization Patterns

| # | Pattern | Location | Description |
|---|---------|----------|-------------|
| 1 | Bit-reversal permutation | `four1.hpp` | Direct bit-reverse index computation per element (`FOR_ALL`), guard `rev > ci` ensures each swap executes on exactly one thread. |
| 2 | Butterfly stages | `four1.hpp` | Twiddle factors precomputed per stage into `CArrayKokkos`, then all butterflies dispatched via `FOR_ALL`. `MATAR_FENCE()` between stages. |
| 3 | Pack / unpack | `twofft.hpp` | Pack is embarrassingly parallel (`DO_ALL`). Unpack uses conjugate symmetry — each thread writes to positions `j` and `n+2-j`, non-overlapping across threads. |
| 4 | Recombination | `realft.hpp` | Twiddle factors precomputed, then `FOR_ALL` over `n/2` independent conjugate pairs. |
| 5 | Pointwise multiply | `convlv.hpp` | Complex multiply/divide in frequency domain — embarrassingly parallel (`DO_ALL`). Also parallel: wrap-around copy and zero-fill of response. |

## Data Representation

Fortran's interleaved `COMPLEX` layout is preserved: a complex array of `n`
elements is stored as `2*n` doubles where element `k` has real part at index
`2k-1` and imaginary part at `2k` (1-based `DFMatrixKokkos<double>`).  This
avoids any layout conversion between routines.

## Files

```
four1.hpp       Parallel Cooley-Tukey radix-2 FFT
realft.hpp      Parallel real-sequence FFT (wraps four1)
twofft.hpp      Parallel simultaneous FFT of two real sequences (wraps four1)
convlv.hpp      Parallel convolution via FFT (wraps twofft + realft)
main.cpp        Driver — signal setup, timing, FFT vs direct comparison
CMakeLists.txt  FetchContent build pulling Kokkos 4.5.01 + MATAR
```

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j$(nproc)
```

Backend options (add to the `cmake` line):

| Option | Effect |
|--------|--------|
| `-DENABLE_OPENMP=ON` | Kokkos OpenMP backend |
| `-DENABLE_CUDA=ON` | Kokkos CUDA backend |
| `-DENABLE_HIP=ON` | Kokkos HIP backend |

## Run

```bash
./convlv        # default N=16 (matches Fortran demo)
./convlv 1024   # larger problem size for benchmarking
```

`N` must be a power of 2 and ≥ 16.

## Expected Output (N=16)

```
   I      CONVLV      Expected
   1       0.000000     0.000000
   2       1.000000     1.000000
   3       1.000000     1.000000
   4       1.000000     1.000000
   5       1.000000     1.000000
   6       1.000000     1.000000
   7      -0.000000     0.000000
   8       1.000000     1.000000
   9       2.000000     2.000000
  10       3.000000     3.000000
  11       3.000000     3.000000
  12       3.000000     3.000000
  13       2.000000     2.000000
  14       1.000000     1.000000
  15       0.000000     0.000000
  16       0.000000     0.000000

  Max |FFT - direct| = 3.30e-15
```

The CONVLV and Expected columns match to at least 6 decimal places (the
Fortran demo's display precision).  The max error is near double-precision
machine epsilon.

## Fortran Reference

```bash
cd fortran/13_spectral_analysis/convlv
make && ./convlv
```
