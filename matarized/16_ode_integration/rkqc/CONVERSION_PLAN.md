---
example: rkqc
chapter: 16_ode_integration
chapter_title: "ODE Integration"
status: not_started
complexity: high
conversion_order: 36
priority: 3
tags: [ode, differential-equation, time-integration, cross-chapter]
dependencies: [bessj, bessj0, bessj1, rk4]
reverse_dependencies: [odeint, shoot, shootf]
---

# RKQC -- ODE Integration

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `rkqc` |
| **Chapter** | 16 -- ODE Integration |
| **Purpose** | Adaptive step-size fifth-order Runge-Kutta-Cash-Karp ODE step. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 44 |
| **Subroutine** | `RKQC` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/16_ode_integration/rkqc/rkqc.f` (44 lines)
- **Driver/demo:** `fortran/16_ode_integration/rkqc/rkqc.dem`
- **Target:** `matarized/16_ode_integration/rkqc/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bessj` (06_special_functions)
  - `bessj0` (06_special_functions)
  - `bessj1` (06_special_functions)
  - `rk4` (16_ode_integration)

### Diagram

```mermaid
graph TD
    rkqc --> bessj
    rkqc --> bessj0
    rkqc --> bessj1
    rkqc --> rk4
```

### Cross-Chapter Dependencies

- `bessj` from chapter 06
- `bessj0` from chapter 06
- `bessj1` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  - `odeint` (16_ode_integration)
  - `shoot` (17_boundary_value_problems)
  - `shootf` (17_boundary_value_problems)

> **Conversion note:** This routine is depended on by 3 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `DERIVS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `DYDX` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `DYSAV` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `EPS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ERRCON` | `REAL` | (scalar) | constant | `constexpr double ERRCON = 6.E-4;` | constant = 6.E-4 |
| `FCOR` | `REAL` | (scalar) | constant | `constexpr double FCOR = .0666666667;` | constant = .0666666667 |
| `HDID` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HNEXT` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HTRY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 10;` | constant = 10 |
| `ONE` | `REAL` | (scalar) | constant | `constexpr double ONE = 1.;` | constant = 1. |
| `SAFETY` | `REAL` | (scalar) | constant | `constexpr double SAFETY = 0.9;` | constant = 0.9 |
| `X` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `YSAV` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `YSCAL` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `YTEMP` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 12  I=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K3: DO 13  I=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None


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
inline void rkqc(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx, int n, double x, double htry, double eps, DFMatrixKokkos<double>& yscal, double hdid, double hnext, double derivs)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `RKQC` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(rkqc_matar_parallel CXX)

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

# --- Build the RKQC example ---
add_executable(rkqc main.cpp)
target_link_libraries(rkqc matar_lib)
target_include_directories(rkqc PRIVATE
    ${SPECIALFUNCTIONS_DIR}/bessj
    ${SPECIALFUNCTIONS_DIR}/bessj0
    ${SPECIALFUNCTIONS_DIR}/bessj1
    ${ODEINTEGRATION_DIR}/rk4
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.
- **Hierarchical parallelism:** For deeply nested loops, consider `FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/16_ode_integration/rkqc
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/16_ode_integration/rkqc
mkdir -p build && cd build
cmake .. && make
./rkqc > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/16_ode_integration/rkqc/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/16_ode_integration/rkqc
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./rkqc > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./rkqc > omp4_output.txt 2>&1
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
| **Conversion order** | 36 of 202 |
| **Priority score** | 3 (reverse dependency count) |
| **Estimated effort** | high (44 Fortran LOC, 4 dependencies) |
| **Prerequisite conversions** | `bessj`, `bessj0`, `bessj1`, `rk4` |
| **Tags** | `ode`, `differential-equation`, `time-integration`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
