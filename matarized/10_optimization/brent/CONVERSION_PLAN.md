---
example: brent
chapter: 10_optimization
chapter_title: "Optimization"
status: not_started
complexity: high
conversion_order: 27
priority: 4
tags: [optimization, minimization, cross-chapter]
dependencies: [bessj0, bessj1, mnbrak]
reverse_dependencies: [dfpmin, frprmn, linmin, powell]
---

# BRENT -- Optimization

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `brent` |
| **Chapter** | 10 -- Optimization |
| **Purpose** | Minimize a function of one variable by Brent's method (parabolic interpolation). |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 79 |
| **Subroutine** | `BRENT` (function) |

## 2. Source Files

- **Fortran source:** `fortran/10_optimization/brent/brent.f` (79 lines)
- **Driver/demo:** `fortran/10_optimization/brent/brent.dem`
- **Target:** `matarized/10_optimization/brent/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bessj0` (06_special_functions)
  - `bessj1` (06_special_functions)
  - `mnbrak` (10_optimization)

### Diagram

```mermaid
graph TD
    brent --> bessj0
    brent --> bessj1
    brent --> mnbrak
```

### Cross-Chapter Dependencies

- `bessj0` from chapter 06
- `bessj1` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  - `dfpmin` (10_optimization)
  - `frprmn` (10_optimization)
  - `linmin` (10_optimization)
  - `powell` (10_optimization)

> **Conversion note:** This routine is depended on by 4 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `AX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `BX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `CGOLD` | `REAL` | (scalar) | constant | `constexpr double CGOLD = .3819660;` | constant = .3819660 |
| `CX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `F` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ITMAX` | `INTEGER` | (scalar) | constant | `constexpr int ITMAX = 100;` | constant = 100 |
| `TOL` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `XMIN` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ZEPS` | `REAL` | (scalar) | constant | `constexpr double ZEPS = 1.0E-10;` | constant = 1.0E-10 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  ITER=1,ITMAX

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
inline double brent(double ax, double bx, double cx, double f, double tol, double xmin)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `BRENT` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(brent_matar_parallel CXX)

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

# --- Build the BRENT example ---
add_executable(brent main.cpp)
target_link_libraries(brent matar_lib)
target_include_directories(brent PRIVATE
    ${SPECIALFUNCTIONS_DIR}/bessj0
    ${SPECIALFUNCTIONS_DIR}/bessj1
    ${OPTIMIZATION_DIR}/mnbrak
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/10_optimization/brent
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/10_optimization/brent
mkdir -p build && cd build
cmake .. && make
./brent > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/10_optimization/brent/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/10_optimization/brent
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./brent > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./brent > omp4_output.txt 2>&1
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
| **Conversion order** | 27 of 202 |
| **Priority score** | 4 (reverse dependency count) |
| **Estimated effort** | high (79 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `bessj0`, `bessj1`, `mnbrak` |
| **Tags** | `optimization`, `minimization`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
