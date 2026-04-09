# Numerical Recipes: Fortran 77 to Performance-Portable C++ with MATAR

This repository stores the classic Numerical Recipes examples in Fortran 77 alongside their performance-portable C++ equivalents built with [MATAR](https://github.com/lanl/MATAR). It serves as a testbed for creating an Agentic AI system that converts legacy Fortran code into modern, performance-portable C++ using MATAR's data structures and parallel dispatch macros.

## Repository Structure

```
nrf77/
├── fortran/          # Original Fortran 77 sources (.f) and demos (.dem)
│   ├── 01_dates_and_calendars/
│   ├── 02_linear_algebra/
│   ├── 03_interpolation/
│   ├── 04_integration/
│   ├── 05_evaluation_of_functions/
│   ├── 06_special_functions/
│   ├── 07_random_numbers/
│   ├── 08_sorting/
│   ├── 09_root_finding/
│   ├── 10_optimization/
│   ├── 11_eigensystems/
│   ├── 12_fourier_transform/
│   ├── 13_spectral_analysis/
│   ├── 14_statistics/
│   ├── 15_curve_fitting/
│   ├── 16_ode_integration/
│   ├── 17_boundary_value_problems/
│   ├── 19_partial_differential_equations/
│   └── data/                              # Shared .dat files used by demos
├── matarized/        # Performance-portable C++ equivalents using MATAR
│   ├── 01_dates_and_calendars/
│   ├── 02_linear_algebra/
│   ├── ...                                # Mirrors fortran/ layout
│   └── 19_partial_differential_equations/
├── MATAR_LLM_CONTEXT.md  # Comprehensive MATAR API and conversion reference
├── readme.md
├── readme.doc             # Original Numerical Recipes diskette documentation
└── names.doc              # Original file listing
```

## Chapter Index

| Dir | Topic | Programs |
| --- | ----- | -------- |
| `01_dates_and_calendars` | Calendar computations | flmoon, julday, badluk, caldat |
| `02_linear_algebra` | Linear equation solution | gaussj, ludcmp, lubksb, tridag, mprove, svbksb, svdcmp, vander, toeplz, sparse |
| `03_interpolation` | Interpolation and extrapolation | polint, ratint, spline, splint, locate, hunt, polcoe, polcof, polin2, bcucof, bcuint, splie2, splin2 |
| `04_integration` | Numerical integration | trapzd, qtrap, qsimp, qromb, midpnt, qromo, midinf, qgaus, gauleg, quad3d |
| `05_evaluation_of_functions` | Series, polynomials, Chebyshev | eulsum, ddpoly, poldiv, chebft, chebev, chder, chint, chebpc, pcshft |
| `06_special_functions` | Gamma, Bessel, error functions, etc. | gammln, factrl, bico, factln, beta, gammp, gammq, gser, gcf, erf, erfc, erfcc, betai, betacf, bessj0, bessy0, bessj1, bessy1, bessy, bessj, bessi0, bessk0, bessi1, bessk1, bessk, bessi, plgndr, sncndn, cel, el2 |
| `07_random_numbers` | Random number generation | ran0, ran1, ran2, ran3, expdev, gasdev, gamdev, poidev, bnldev, irbit1, irbit2, ran4, des, desks |
| `08_sorting` | Sorting and ranking | piksrt, piksr2, shell, sort, sort2, indexx, sort3, rank, eclass, eclazz, qcksrt, mdian1, mdian2 |
| `09_root_finding` | Root finding and nonlinear equations | scrsho, zbrac, zbrak, rtbis, rtflsp, rtsec, zbrent, rtnewt, rtsafe, laguer, zroots, qroot, mnewt |
| `10_optimization` | Minimization and maximization | mnbrak, golden, brent, dbrent, amoeba, powell, linmin, f1dim, frprmn, df1dim, dfpmin, simplx, simp1, simp2, simp3, anneal, link |
| `11_eigensystems` | Eigenvalue problems | jacobi, eigsrt, tred2, tqli, balanc, elmhes, hqr |
| `12_fourier_transform` | Fast Fourier Transform | four1, twofft, realft, sinft, cosft, fourn |
| `13_spectral_analysis` | Fourier and spectral applications | convlv, correl, spctrm, memcof, fixrts, predic, evlmem, smooft |
| `14_statistics` | Statistical tests and descriptors | moment, ttest, avevar, tutest, tptest, ftest, chsone, chstwo, ksone, kstwo, probks, cntab1, cntab2, pearsn, spear, crank, kendl1, kendl2 |
| `15_curve_fitting` | Least-squares and robust fitting | fit, lfit, covsrt, svdfit, svdvar, fpoly, fleg, mrqmin, mrqcof, fgauss, medfit, rofunc |
| `16_ode_integration` | Ordinary differential equations | rk4, rkdumb, odeint, mmid, bsstep, pzextr, rzextr, rkqc |
| `17_boundary_value_problems` | Two-point boundary value problems | shoot, shootf, solvde, bksub, pinvs, red, sfroid, difeq |
| `19_partial_differential_equations` | Elliptic PDEs | sor, adi |

## About MATAR

[MATAR](https://github.com/lanl/MATAR) (MATrix and ARray) is a header-only C++ library that provides performance-portable multi-dimensional data structures, parallel loop macros, and reduction operations. MATAR supports CPU serial, OpenMP, Pthreads, CUDA, and HIP backends from a single source code — the developer writes to MATAR's API and the backend is selected at build time. Under the hood MATAR can use Kokkos for device execution, but user code interacts exclusively with MATAR types and macros.

### Data Structure Taxonomy

MATAR types are named along four axes:

| Axis | Options | Meaning |
| ---- | ------- | ------- |
| **Layout** | `C` / `F` | Row-major (C-style) or column-major (Fortran-style) |
| **Index base** | `Array` / `Matrix` | 0-based `[0, N)` or 1-based `[1, N]` |
| **Residence** | `Host` / `Device` / `Dual` | CPU-only, device-only, or both with explicit sync |
| **Ownership** | (none) / `View` | Allocates memory or wraps an existing pointer |

These combine into names like `CArrayDevice`, `FMatrixDual`, `ViewCArrayHost`, etc. For Fortran interoperability the key types are `FMatrixDevice` and `ViewFMatrixDual` — column-major, 1-based, matching Fortran's native layout.

### Parallel Dispatch

| Macro | Range Style | Purpose |
| ----- | ----------- | ------- |
| `FOR_ALL` | Half-open `[lo, hi)` | Parallel loop (1D/2D/3D) |
| `DO_ALL` | Inclusive `[lo, hi]` | Parallel loop for 1-based indexing |
| `FOR_REDUCE_SUM/MAX/MIN/PRODUCT` | Half-open | Parallel reductions |
| `DO_REDUCE_SUM/MAX/MIN` | Inclusive | Parallel reductions for 1-based indexing |
| `RUN` | Single iteration | Execute once on device |

All MATAR programs require `MATAR_INITIALIZE` / `MATAR_FINALIZE` bracketing, and `MATAR_FENCE()` between kernels that have a data dependency.

### Comprehensive Reference

See [`MATAR_LLM_CONTEXT.md`](MATAR_LLM_CONTEXT.md) for the full API reference, including detailed construction patterns, sparse/ragged types, hierarchical parallelism, host/device transfer rules, fence placement, MPI communication plans, and device kernel constraints.

## Conversion Approach

Each Fortran subroutine in `fortran/<chapter>/<name>.f` has a corresponding C++ file at `matarized/<chapter>/<name>.cpp`. The conversion follows a systematic process:

### 1. Data Structure Mapping

| Fortran Pattern | MATAR Replacement | Notes |
| --------------- | ----------------- | ----- |
| `REAL*8 A(N,M)` | `FMatrixDevice<double>(N, M)` | Column-major, 1-based — preserves Fortran semantics |
| `INTEGER A(N)` | `FMatrixDevice<int>(N)` | Or `CArrayDevice<int>(N)` when switching to 0-based |
| Subroutine array argument | `ViewFMatrixDual<double>(ptr, N, M)` | Non-owning view wrapping the caller's memory |
| `COMMON` block arrays | Struct of MATAR arrays or `Dual` types | Shared state across translation units |

### 2. Loop Mapping

| Fortran | MATAR | Notes |
| ------- | ----- | ----- |
| `DO i = 1, N` ... `ENDDO` | `DO_ALL(i, 1, N, { ... });` | Inclusive range, 1-based |
| Nested `DO` with accumulator | `DO_REDUCE_SUM(i, 1, N, loc, { loc += ...; }, total);` | Parallel reduction |
| `DO` with data dependency on inner dim | `DO_ALL` outer + plain `for` inner | Keep serial dimensions serial |

### 3. Correctness and Synchronization

- `MATAR_FENCE()` between kernels where the second reads data written by the first
- `update_device()` after host writes, before device reads
- `update_host()` after device writes, before host reads
- Reduction results are available immediately (implicit fence)

### 4. Algorithmic Preservation

The conversions maintain the same algorithmic structure as the original Fortran for traceability. Numerical output is validated against the original Fortran `.dem` driver programs.

## Build Configuration

### CMake

```cmake
find_package(Matar REQUIRED)
target_link_libraries(my_target matar)
```

### Backend Selection

| Define | Backend |
| ------ | ------- |
| `HAVE_KOKKOS` | Enables MATAR's device/dual types and parallel macros (via Kokkos backend) |
| `HAVE_CUDA` | CUDA GPU execution |
| `HAVE_HIP` | AMD GPU execution |
| `HAVE_OPENMP` | OpenMP threading |
| `HAVE_THREADS` | Pthreads |

Without a GPU backend, MATAR defaults to CPU execution (serial or OpenMP depending on build flags).

### Running

```bash
# Serial
./my_program

# OpenMP
export OMP_NUM_THREADS=8 && ./my_program

# Threaded (HAVE_THREADS)
./my_program 
```

## Original Sources

The Fortran 77 programs are from *Numerical Recipes: The Art of Scientific Computing* by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling, published by Cambridge University Press. The `.dem` files are demonstration/driver programs; the `.f` files are the subroutine implementations.
