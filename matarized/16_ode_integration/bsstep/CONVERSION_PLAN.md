---
example: bsstep
chapter: 16_ode_integration
chapter_title: "ODE Integration"
status: not_started
complexity: high
conversion_order: 182
priority: 0
tags: [ode, differential-equation, time-integration, cross-chapter]
dependencies: [bessj, bessj0, bessj1, mmid, rzextr]
reverse_dependencies: []
---

# BSSTEP -- ODE Integration

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `bsstep` |
| **Chapter** | 16 -- ODE Integration |
| **Purpose** | Bulirsch-Stoer ODE step with Richardson extrapolation. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 39 |
| **Subroutine** | `BSSTEP` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/16_ode_integration/bsstep/bsstep.f` (39 lines)
- **Driver/demo:** `fortran/16_ode_integration/bsstep/bsstep.dem`
- **Target:** `matarized/16_ode_integration/bsstep/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bessj` (06_special_functions)
  - `bessj0` (06_special_functions)
  - `bessj1` (06_special_functions)
  - `mmid` (16_ode_integration)
  - `rzextr` (16_ode_integration)

### Diagram

```mermaid
graph TD
    bsstep --> bessj
    bsstep --> bessj0
    bsstep --> bessj1
    bsstep --> mmid
    bsstep --> rzextr
```

### Cross-Chapter Dependencies

- `bessj` from chapter 06
- `bessj0` from chapter 06
- `bessj1` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `DERIVS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `DYDX` | `REAL` | NV | parameter (input) | `DFMatrixKokkos<double>(NV)` |  |
| `DYSAV` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `EPS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `GROW` | `REAL` | (scalar) | constant | `constexpr double GROW = 1.2E0;` | constant = 1.2E0 |
| `HDID` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HNEXT` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HTRY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `IMAX` | `INTEGER` | (scalar) | constant | `constexpr int IMAX = 11;` | constant = 11 |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 10;` | constant = 10 |
| `NSEQ` | `INTEGER` | IMAX | local | `DFMatrixKokkos<int>(IMAX)` |  |
| `NUSE` | `INTEGER` | (scalar) | constant | `constexpr int NUSE = 7;` | constant = 7 |
| `NV` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `ONE` | `REAL` | (scalar) | constant | `constexpr double ONE = 1.E0;` | constant = 1.E0 |
| `SHRINK` | `REAL` | (scalar) | constant | `constexpr double SHRINK = .95E0;` | constant = .95E0 |
| `X` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL` | NV | parameter (input) | `DFMatrixKokkos<double>(NV)` |  |
| `YERR` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `YSAV` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `YSCAL` | `REAL` | NV | parameter (input) | `DFMatrixKokkos<double>(NV)` |  |
| `YSEQ` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=1,NV

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: X

### K2: DO 10  I=1,IMAX

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: X

### K3: DO 12  J=1,NV

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: X


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
inline void bsstep(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx, int nv, double x, double htry, double eps, DFMatrixKokkos<double>& yscal, double hdid, double hnext, double derivs)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `BSSTEP` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(bsstep_matar_parallel CXX)

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

# --- Build the BSSTEP example ---
add_executable(bsstep main.cpp)
target_link_libraries(bsstep matar_lib)
target_include_directories(bsstep PRIVATE
    ${SPECIALFUNCTIONS_DIR}/bessj
    ${SPECIALFUNCTIONS_DIR}/bessj0
    ${SPECIALFUNCTIONS_DIR}/bessj1
    ${ODEINTEGRATION_DIR}/mmid
    ${ODEINTEGRATION_DIR}/rzextr
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
cd fortran/16_ode_integration/bsstep
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/16_ode_integration/bsstep
mkdir -p build && cd build
cmake .. && make
./bsstep > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/16_ode_integration/bsstep/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/16_ode_integration/bsstep
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./bsstep > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./bsstep > omp4_output.txt 2>&1
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
| **Conversion order** | 182 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (39 Fortran LOC, 5 dependencies) |
| **Prerequisite conversions** | `bessj`, `bessj0`, `bessj1`, `mmid`, `rzextr` |
| **Tags** | `ode`, `differential-equation`, `time-integration`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
