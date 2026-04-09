# SOR -- Successive Over-Relaxation (Fortran to MATAR Conversion)

This directory contains a three-step conversion of the Numerical Recipes
Fortran 77 SOR solver (`sor.f` + `sor.dem`) to performance-portable C++
using MATAR.

## Algorithm

SOR solves a 2D elliptic PDE on a JMAX x JMAX grid using red-black
(checkerboard) ordering.  In each iteration, only cells of one color are
updated; their stencil neighbors are all of the opposite color, making all
updates within a color independent.  This property is exploited in step 2
for full parallelism.  Chebyshev acceleration of the relaxation parameter
`omega` speeds convergence.

## Conversion Steps

| Step | Directory | Data Types | Parallelism | Dependencies |
|------|-----------|-----------|-------------|-------------|
| 0 | `step0_cpp_baseline/` | Raw C++ arrays (0-based) | None | None |
| 1 | `step1_matar_serial/` | `FMatrix<double>` (1-based, host) | None | MATAR (header-only) |
| 2 | `step2_matar_parallel/` | `DFMatrixKokkos<double>` (1-based, dual) | `DO_REDUCE_SUM` | MATAR + Kokkos |

Each step is self-contained with its own `CMakeLists.txt`.  All three
produce identical output.

## Building

### Step 0 -- C++ Baseline

```bash
cd step0_cpp_baseline
cmake -B build
cmake --build build
./build/sor
```

### Step 1 -- MATAR Serial (host only, no Kokkos)

```bash
cd step1_matar_serial
cmake -B build
cmake --build build
./build/sor
```

### Step 2 -- MATAR Parallel (Kokkos backend)

```bash
cd step2_matar_parallel

# Serial backend (default)
cmake -B build
cmake --build build
./build/sor

# OpenMP backend
cmake -B build-omp -DENABLE_OPENMP=ON
cmake --build build-omp
./build-omp/sor

# CUDA backend
cmake -B build-cuda -DENABLE_CUDA=ON
cmake --build build-cuda
./build-cuda/sor
```

## Expected Output

All steps produce the same output (matching the original Fortran demo):

```
 SOR Solution:
   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
   0.00 -0.02 -0.04 -0.06 -0.08 -0.09 -0.08 -0.06 -0.04 -0.02  0.00
   0.00 -0.04 -0.09 -0.13 -0.17 -0.19 -0.17 -0.13 -0.09 -0.04  0.00
   0.00 -0.06 -0.13 -0.20 -0.28 -0.32 -0.28 -0.20 -0.13 -0.06  0.00
   0.00 -0.08 -0.17 -0.28 -0.41 -0.55 -0.41 -0.28 -0.17 -0.08  0.00
   0.00 -0.09 -0.19 -0.32 -0.55 -1.05 -0.55 -0.32 -0.19 -0.09  0.00
   0.00 -0.08 -0.17 -0.28 -0.41 -0.55 -0.41 -0.28 -0.17 -0.08  0.00
   0.00 -0.06 -0.13 -0.20 -0.28 -0.32 -0.28 -0.20 -0.13 -0.06  0.00
   0.00 -0.04 -0.09 -0.13 -0.17 -0.19 -0.17 -0.13 -0.09 -0.04  0.00
   0.00 -0.02 -0.04 -0.06 -0.08 -0.09 -0.08 -0.06 -0.04 -0.02  0.00
   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00

 Test that solution satisfies Difference Eqns:
         0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
         0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00
         0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00  0.00
         0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00
         0.00  0.00 -0.00  0.00  2.00  0.00 -0.00  0.00  0.00
         0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00
         0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00  0.00
         0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00 -0.00  0.00
         0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
```

The verification grid should be zero everywhere except at the center point
`(6,6)` where the source `F(6,6) = 2.0` was applied.

## Parallelization Notes

The red-black coloring in SOR naturally maps to MATAR's parallel dispatch:

- `DO_REDUCE_SUM` parallelizes the 2D interior grid `[2, jmax-1] x [2, jmax-1]`
- A parity check `(j + l) % 2 == n % 2` selects the active color
- Within each color, writes to `u(j, l)` are conflict-free
- The residual norm `anorm` is accumulated via the reduction variable

For larger grids, this parallelism maps directly to GPU threads via the
Kokkos CUDA or HIP backends with no code changes.
