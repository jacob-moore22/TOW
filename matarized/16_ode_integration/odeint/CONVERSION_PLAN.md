---
example: odeint
chapter: 16_ode_integration
chapter_title: "ODE Integration"
status: not_started
complexity: high
conversion_order: 56
priority: 2
tags: [ode, differential-equation, time-integration, cross-chapter]
dependencies: [bessj, bessj0, bessj1, rk4, rkqc]
reverse_dependencies: [shoot, shootf]
---

# ODEINT -- ODE Integration

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `odeint` |
| **Chapter** | 16 -- ODE Integration |
| **Purpose** | Driver for ODE integration with adaptive step-size control. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 58 |
| **Subroutine** | `ODEINT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/16_ode_integration/odeint/odeint.f` (58 lines)
- **Driver/demo:** `fortran/16_ode_integration/odeint/odeint.dem`
- **Target:** `matarized/16_ode_integration/odeint/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bessj` (06_special_functions)
  - `bessj0` (06_special_functions)
  - `bessj1` (06_special_functions)
  - `rk4` (16_ode_integration)
  - `rkqc` (16_ode_integration)

### Diagram

```mermaid
graph TD
    odeint --> bessj
    odeint --> bessj0
    odeint --> bessj1
    odeint --> rk4
    odeint --> rkqc
```

### Cross-Chapter Dependencies

- `bessj` from chapter 06
- `bessj0` from chapter 06
- `bessj1` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  - `shoot` (17_boundary_value_problems)
  - `shootf` (17_boundary_value_problems)

> **Conversion note:** This routine is depended on by 2 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `DERIVS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `DYDX` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `EPS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `H1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HMIN` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `MAXSTP` | `INTEGER` | (scalar) | constant | `constexpr int MAXSTP = 10000;` | constant = 10000 |
| `NBAD` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 10;` | constant = 10 |
| `NOK` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NVAR` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `RK QC` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `TINY` | `REAL` | (scalar) | constant | `constexpr double TINY = 1.E-30;` | constant = 1.E-30 |
| `TWO` | `REAL` | (scalar) | constant | `constexpr double TWO = 2.0;` | constant = 2.0 |
| `X1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `YSCAL` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `YSTART` | `REAL` | NVAR | parameter (input) | `DFMatrixKokkos<double>(NVAR)` |  |
| `ZERO` | `REAL` | (scalar) | constant | `constexpr double ZERO = 0.0;` | constant = 0.0 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=1,NVAR

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.

### K2: DO 16  NSTP=1,MAXSTP

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.

### K3: DO 12  I=1,NVAR

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.

### K4: DO 13  I=1,NVAR

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.

### K5: DO 14  I=1,NVAR

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.

### K6: DO 15  I=1,NVAR

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: KOUNT, NBAD, NOK Array write(s) not indexed by loop variable: XP. Verify thread safety.


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
inline void odeint(DFMatrixKokkos<double>& ystart, int nvar, double x1, double x2, double eps, double h1, double hmin, int nok, int nbad, double derivs, double rk qc)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `ODEINT` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(odeint_matar_parallel CXX)

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
set(SPECIALFUNCTIONS_DIR ${MATARIZED_ROOT}/06_special_functions)
set(ODEINTEGRATION_DIR   ${MATARIZED_ROOT}/16_ode_integration)

# --- Build the ODEINT example ---
add_executable(odeint main.cpp)
target_link_libraries(odeint matar_lib)
target_include_directories(odeint PRIVATE
    ${SPECIALFUNCTIONS_DIR}/bessj
    ${SPECIALFUNCTIONS_DIR}/bessj0
    ${SPECIALFUNCTIONS_DIR}/bessj1
    ${ODEINTEGRATION_DIR}/rk4
    ${ODEINTEGRATION_DIR}/rkqc
)
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
cd fortran/16_ode_integration/odeint
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/16_ode_integration/odeint
mkdir -p build && cd build
cmake .. && make
./odeint > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/16_ode_integration/odeint/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/16_ode_integration/odeint
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./odeint > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./odeint > omp4_output.txt 2>&1
# Verify: omp1 output must exactly match serial output
diff serial_output.txt omp1_output.txt
# Verify: omp4 output must match within floating-point tolerance
```


### Pass Criteria

- Max absolute difference vs. Fortran reference: **< 1e-10** (double precision)

- OpenMP results must be deterministic across repeated runs

- No runtime errors, memory leaks, or Kokkos warnings


## 11. Agent Metadata

| Field | Value |
|-------|-------|
| **Conversion order** | 56 of 202 |
| **Priority score** | 2 (reverse dependency count) |
| **Estimated effort** | high (58 Fortran LOC, 5 dependencies) |
| **Prerequisite conversions** | `bessj`, `bessj0`, `bessj1`, `rk4`, `rkqc` |
| **Tags** | `ode`, `differential-equation`, `time-integration`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
