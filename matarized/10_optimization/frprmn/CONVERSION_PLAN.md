---
example: frprmn
chapter: 10_optimization
chapter_title: "Optimization"
status: not_started
complexity: high
conversion_order: 190
priority: 0
tags: [optimization, minimization, cross-chapter]
dependencies: [bessj0, bessj1, brent, f1dim, linmin, mnbrak]
reverse_dependencies: []
---

# FRPRMN -- Optimization

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `frprmn` |
| **Chapter** | 10 -- Optimization |
| **Purpose** | Minimize in N-D using Fletcher-Reeves-Polak-Ribiere conjugate gradient. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 35 |
| **Subroutine** | `FRPRMN` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/10_optimization/frprmn/frprmn.f` (35 lines)
- **Driver/demo:** `fortran/10_optimization/frprmn/frprmn.dem`
- **Target:** `matarized/10_optimization/frprmn/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bessj0` (06_special_functions)
  - `bessj1` (06_special_functions)
  - `brent` (10_optimization)
  - `f1dim` (10_optimization)
  - `linmin` (10_optimization)
  - `mnbrak` (10_optimization)

### Diagram

```mermaid
graph TD
    frprmn --> bessj0
    frprmn --> bessj1
    frprmn --> brent
    frprmn --> f1dim
    frprmn --> linmin
    frprmn --> mnbrak
```

### Cross-Chapter Dependencies

- `bessj0` from chapter 06
- `bessj1` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `EPS` | `REAL` | (scalar) | constant | `constexpr double EPS = 1.E-10;` | constant = 1.E-10 |
| `FRET` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `FTOL` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `G` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `H` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `ITER` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `ITMAX` | `INTEGER` | (scalar) | constant | `constexpr int ITMAX = 200;` | constant = 200 |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 50;` | constant = 50 |
| `P` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `XI` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: DGG, GG

### K2: DO 14  ITS=1,ITMAX

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: DGG, GG

### K3: DO 12  J=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: DGG, GG

### K4: DO 13  J=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: DGG, GG


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
inline void frprmn(DFMatrixKokkos<double>& p, int n, double ftol, int iter, double fret)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `FRPRMN` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(frprmn_matar_parallel CXX)

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
set(OPTIMIZATION_DIR     ${MATARIZED_ROOT}/10_optimization)

# --- Build the FRPRMN example ---
add_executable(frprmn main.cpp)
target_link_libraries(frprmn matar_lib)
target_include_directories(frprmn PRIVATE
    ${SPECIALFUNCTIONS_DIR}/bessj0
    ${SPECIALFUNCTIONS_DIR}/bessj1
    ${OPTIMIZATION_DIR}/brent
    ${OPTIMIZATION_DIR}/f1dim
    ${OPTIMIZATION_DIR}/linmin
    ${OPTIMIZATION_DIR}/mnbrak
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
cd fortran/10_optimization/frprmn
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/10_optimization/frprmn
mkdir -p build && cd build
cmake .. && make
./frprmn > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/10_optimization/frprmn/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/10_optimization/frprmn
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./frprmn > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./frprmn > omp4_output.txt 2>&1
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
| **Conversion order** | 190 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (35 Fortran LOC, 6 dependencies) |
| **Prerequisite conversions** | `bessj0`, `bessj1`, `brent`, `f1dim`, `linmin`, `mnbrak` |
| **Tags** | `optimization`, `minimization`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
