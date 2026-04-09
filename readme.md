# Numerical Recipes: Fortran 77 to Performance-Portable C++ with MATAR

This repository stores the classic Numerical Recipes examples in Fortran 77 alongside their performance-portable C++ equivalents built with [MATAR](https://github.com/lanl/MATAR). It serves as a testbed for creating an Agentic AI system that converts legacy Fortran code into modern, performance-portable C++ using Kokkos-backed data structures.

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
├── readme.md
├── readme.doc        # Original Numerical Recipes diskette documentation
└── names.doc         # Original file listing
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

[MATAR](https://github.com/lanl/MATAR) is a C++ library for simple, fast, and memory-efficient multi-dimensional data representations that are portable across CPUs and GPUs via Kokkos. Each C++ file under `matarized/` will use MATAR's data types (CArray, FMatrix, ViewCArray, etc.) to replace Fortran arrays and achieve performance portability.

## Conversion Approach

Each Fortran subroutine in `fortran/<chapter>/<name>.f` has a corresponding C++ file at `matarized/<chapter>/<name>.cpp`. The conversion targets:

1. Replace Fortran arrays with MATAR types (CArray, FMatrix, etc.)
2. Replace DO loops with Kokkos parallel dispatch where beneficial
3. Preserve numerical correctness against the original Fortran output
4. Maintain the same algorithmic structure for traceability

## Original Sources

The Fortran 77 programs are from *Numerical Recipes: The Art of Scientific Computing* by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling, published by Cambridge University Press. The `.dem` files are demonstration/driver programs; the `.f` files are the subroutine implementations.
