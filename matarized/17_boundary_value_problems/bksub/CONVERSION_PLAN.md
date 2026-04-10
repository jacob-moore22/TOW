---
example: bksub
chapter: 17_boundary_value_problems
chapter_title: "Boundary Value Problems"
status: library_only
complexity: low
conversion_order: 58
priority: 1
tags: [boundary-value, differential-equation, library, leaf]
dependencies: []
reverse_dependencies: [sfroid]
---

# BKSUB -- Boundary Value Problems

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `bksub` |
| **Chapter** | 17 -- Boundary Value Problems |
| **Purpose** | Back-substitution step for the banded system in the relaxation solver. |
| **Status** | `library_only` |
| **Complexity** | `low` |
| **Fortran LOC** | 24 |
| **Subroutine** | `BKSUB` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/17_boundary_value_problems/bksub/bksub.f` (24 lines)
- **Driver/demo:** _(none -- library-only or program-in-.f)_
- **Target:** `matarized/17_boundary_value_problems/bksub/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    bksub
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
| `JF` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `K1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `K2` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NB` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCI` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCJ` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NCK` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NE` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 13  K=K2,K1, step -1

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 12  J=1,NBF

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K3: DO 11  I=1,NE

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K4: DO 16  K=K1,K2

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K5: DO 14  I=1,NB

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K6: DO 15  I=1,NBF

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
inline void bksub(int ne, int nb, int jf, int k1, int k2, DFMatrixKokkos<double>& c, int nci, int ncj, int nck)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `BKSUB` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(bksub_matar_parallel CXX)

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


# --- Build the BKSUB example ---
add_executable(bksub main.cpp)
target_link_libraries(bksub matar_lib)

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
cd matarized/17_boundary_value_problems/bksub
mkdir -p build && cd build
cmake .. && make
./bksub > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/17_boundary_value_problems/bksub/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/17_boundary_value_problems/bksub
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./bksub > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./bksub > omp4_output.txt 2>&1
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
| **Conversion order** | 58 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | low (24 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `boundary-value`, `differential-equation`, `library`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
