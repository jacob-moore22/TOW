---
example: sor
chapter: 19_partial_differential_equations
chapter_title: "Partial Differential Equations"
status: converted
complexity: medium
conversion_order: 121
priority: 0
tags: [pde, differential-equation, stencil, leaf]
dependencies: []
reverse_dependencies: []
---

# SOR -- Partial Differential Equations

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `sor` |
| **Chapter** | 19 -- Partial Differential Equations |
| **Purpose** | Successive Over-Relaxation solver for 2D elliptic PDEs with red-black ordering. |
| **Status** | `converted` |
| **Complexity** | `medium` |
| **Fortran LOC** | 36 |
| **Subroutine** | `SOR` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/19_partial_differential_equations/sor/sor.f` (36 lines)
- **Driver/demo:** `fortran/19_partial_differential_equations/sor/sor.dem`
- **Target:** `matarized/19_partial_differential_equations/sor/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    sor
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `B` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `C` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `D` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `E` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `EPS` | `REAL*8` | (scalar) | constant | `constexpr double EPS = 1.D-5;` | constant = 1.D-5 |
| `F` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `HALF` | `REAL*8` | (scalar) | constant | `constexpr double HALF = .5D0;` | constant = .5D0 |
| `JMAX` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `MAXITS` | `INTEGER` | (scalar) | constant | `constexpr int MAXITS = 1000;` | constant = 1000 |
| `ONE` | `REAL*8` | (scalar) | constant | `constexpr double ONE = 1.D0;` | constant = 1.D0 |
| `QTR` | `REAL*8` | (scalar) | constant | `constexpr double QTR = .25D0;` | constant = .25D0 |
| `RJAC` | `REAL*8` | (scalar) | parameter (input) | `double` |  |
| `U` | `REAL*8` | JMAX, JMAX | parameter (input) | `DFMatrixKokkos<double>(JMAX, JMAX)` |  |
| `ZERO` | `REAL*8` | (scalar) | constant | `constexpr double ZERO = 0.D0;` | constant = 0.D0 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  J=2,JMAX-1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: ANORM, ANORMF

### K2: DO 11  L=2,JMAX-1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: ANORM, ANORMF

### K3: DO 15  N=1,MAXITS  (OUTER ITERATION LOOP)

- **Thread safety:** `inherently_serial`
- **Recommended macro:** _serial `for` loop (or `while` with convergence check)_
- **Notes:** Each SOR iteration depends on the previous one's solution and omega value. Must remain serial. Contains the inner sweep kernels K4/K5.

### K4+K5: DO 14 J / DO 13 L  (RED-BLACK INTERIOR SWEEP)

- **Thread safety:** `safe` (with red-black coloring)
- **Recommended macro:** `DO_ALL` for spatial sweep + `DO_REDUCE_SUM` for ANORM
- **Notes:** The `IF(MOD(J+L,2).EQ.MOD(N,2))` guard implements **red-black (checkerboard) ordering**: within a single color pass, all updated cells have neighbors of the opposite color only. This makes all updates within a color independent. **Strategy:** Split into two `DO_ALL` passes per iteration (red then black, or vice versa). Each pass updates only same-color cells; their stencil reads are from the other color which was updated in the previous pass. The ANORM accumulation is a `DO_REDUCE_SUM` (or `DO_REDUCE_MAX` for max-norm). See the existing `step2_matar_parallel/sor.cpp` for the reference implementation.


### Thread-Safety Legend

| Classification | Meaning | Action |
|---------------|---------|--------|
| `safe` | No write conflicts | Parallelize directly with `DO_ALL` |
| `reduction` | Accumulation to scalar | Use `DO_REDUCE_SUM` / `DO_REDUCE_MAX` |
| `unsafe_review` | Potential race condition | Restructure: inner serial loop or phased approach |
| `inherently_serial` | Sequential data dependency | Keep as serial `for` inside parallel region |

## 7. Conversion Strategy

### Proposed C++ Signature

```cpp
inline void sor(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b, DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& e, DFMatrixKokkos<double>& f, DFMatrixKokkos<double>& u, int jmax, double rjac)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `SOR` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(sor_matar_parallel CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# --- Kokkos backend selection (Serial is always on) ---
set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable Kokkos serial backend")

option(ENABLE_OPENMP "Enable OpenMP backend" OFF)
option(ENABLE_CUDA   "Enable CUDA backend"   OFF)
option(ENABLE_HIP    "Enable HIP backend"    OFF)

if(ENABLE_OPENMP)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
endif()
if(ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA        ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
endif()
if(ENABLE_HIP)
    set(Kokkos_ENABLE_HIP ON CACHE BOOL "")
endif()

# --- Fetch Kokkos ---
FetchContent_Declare(
    kokkos
    GIT_REPOSITORY https://github.com/kokkos/kokkos.git
    GIT_TAG        4.5.01
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(kokkos)

# --- Fetch MATAR (header-only -- bypass its CMakeLists.txt) ---
FetchContent_Declare(
    matar
    GIT_REPOSITORY https://github.com/lanl/MATAR.git
    GIT_TAG        main
    GIT_SHALLOW    TRUE
)
FetchContent_GetProperties(matar)
if(NOT matar_POPULATED)
    FetchContent_Populate(matar)
endif()

add_library(matar_lib INTERFACE)
target_include_directories(matar_lib INTERFACE ${matar_SOURCE_DIR}/src/include)
target_link_libraries(matar_lib INTERFACE Kokkos::kokkos)
target_compile_definitions(matar_lib INTERFACE HAVE_KOKKOS=1)

# --- Cross-chapter dependency headers ---
set(MATARIZED_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)


# --- Build the SOR example ---
add_executable(sor main.cpp)
target_link_libraries(sor matar_lib)

```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Reduction fusion:** If multiple reductions share the same loop bounds, consider fusing them into a single pass to reduce kernel launch overhead.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.
- **Hierarchical parallelism:** For deeply nested loops, consider `FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/19_partial_differential_equations/sor
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/19_partial_differential_equations/sor
mkdir -p build && cd build
cmake .. && make
./sor > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/19_partial_differential_equations/sor/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/19_partial_differential_equations/sor
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./sor > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./sor > omp4_output.txt 2>&1
# Verify: omp1 output must exactly match serial output
diff serial_output.txt omp1_output.txt
# Verify: omp4 output must match within floating-point tolerance
```


### Pass Criteria

- Max absolute difference vs. Fortran reference: **< 1e-10** (double precision)

- OpenMP results must be deterministic across repeated runs

- No runtime errors, memory leaks, or Kokkos warnings


## Conversion Audit (Already Converted)

This example has a complete three-step conversion:

1. **step0_cpp_baseline/** -- Plain C++ translation (128 lines) with `CMakeLists.txt`
2. **step1_matar_serial/** -- MATAR serial using `FMatrix` host types (139 lines)
3. **step2_matar_parallel/** -- Full MATAR parallel with `DFMatrixKokkos` + red-black coloring (167 lines) with `CMakeLists.txt`

### Key Design Decisions

- **Red-black coloring:** The Fortran `MOD(J+L,2).EQ.MOD(N,2)` conditional is restructured into two explicit `DO_ALL` passes per iteration -- one for even-parity cells, one for odd-parity. This eliminates the branch and makes each pass fully independent.
- **Chebyshev acceleration:** The omega update (`OMEGA = 1/(1 - 0.25*RJAC^2*OMEGA)`) remains serial between iterations as it has a sequential dependency.
- **Convergence check:** ANORM uses `DO_REDUCE_SUM` across the 2D interior grid, then compared to `ANORMF * EPS` on the host.
- **Data types:** `DFMatrixKokkos<double>` (Dual, Fortran-layout, 1-based) for all seven coefficient/solution arrays.

### Status: COMPLETE

See also: `matarized/19_partial_differential_equations/sor/README.md` for full documentation.

## 11. Agent Metadata

| Field | Value |
|-------|-------|
| **Conversion order** | 121 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (36 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `pde`, `differential-equation`, `stencil`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
