---
example: pinvs
chapter: 17_boundary_value_problems
chapter_title: "Boundary Value Problems"
status: library_only
complexity: high
conversion_order: 65
priority: 1
tags: [boundary-value, differential-equation, library, leaf]
dependencies: []
reverse_dependencies: [sfroid]
---

# PINVS -- Boundary Value Problems

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `pinvs` |
| **Chapter** | 17 -- Boundary Value Problems |
| **Purpose** | Diagonalize and back-substitute a section of the banded matrix. |
| **Status** | `library_only` |
| **Complexity** | `high` |
| **Fortran LOC** | 61 |
| **Subroutine** | `PINVS` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/17_boundary_value_problems/pinvs/pinvs.f` (61 lines)
- **Driver/demo:** _(none -- library-only or program-in-.f)_
- **Target:** `matarized/17_boundary_value_problems/pinvs/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    pinvs
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `sfroid` (17_boundary_value_problems)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `C` | `REAL` | NCI, NCJ, NCK | parameter (input) | `DFMatrixKokkos<double>(NCI, NCJ, NCK)` |  |
| `IE1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `IE2` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `INDXR` | `INTEGER` | NMAX | local | `DFMatrixKokkos<int>(NMAX)` |  |
| `JC1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `JE1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `JSF` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `K` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCI` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCJ` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCK` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 10;` | constant = 10 |
| `NSI` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NSJ` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `ONE` | `REAL` | (scalar) | constant | `constexpr double ONE = 1.;` | constant = 1. |
| `PSCL` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `S` | `REAL` | NSI, NSJ | parameter (input) | `DFMatrixKokkos<double>(NSI, NSJ)` |  |
| `ZERO` | `REAL` | (scalar) | constant | `constexpr double ZERO = 0.;` | constant = 0. |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  I=IE1,IE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 11  J=JE1,JE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K3: DO 18  ID=IE1,IE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K4: DO 14  I=IE1,IE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K5: DO 13  J=JE1,JE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K6: DO 15  J=JE1,JSF

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K7: DO 17  I=IE1,IE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K8: DO 16  J=JE1,JSF

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K9: DO 21  I=IE1,IE2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K10: DO 19  J=JS1,JSF

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.


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
inline void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k, DFMatrixKokkos<double>& c, int nci, int ncj, int nck, DFMatrixKokkos<double>& s, int nsi, int nsj)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `PINVS` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(pinvs_matar_parallel CXX)

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


# --- Build the PINVS example ---
add_executable(pinvs main.cpp)
target_link_libraries(pinvs matar_lib)

```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.
- **Hierarchical parallelism:** For deeply nested loops, consider `FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.

## 10. Validation Plan

### Reference Output

_No standalone Fortran driver (.dem) exists.  Validate by calling this routine from a dependent example._


### Serial Validation

```bash
cd matarized/17_boundary_value_problems/pinvs
mkdir -p build && cd build
cmake .. && make
./pinvs > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/17_boundary_value_problems/pinvs/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/17_boundary_value_problems/pinvs
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./pinvs > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./pinvs > omp4_output.txt 2>&1
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
| **Conversion order** | 65 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | high (61 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `boundary-value`, `differential-equation`, `library`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
