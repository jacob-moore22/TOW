---
example: link
chapter: 10_optimization
chapter_title: "Optimization"
status: library_only
complexity: high
conversion_order: 63
priority: 1
tags: [optimization, minimization, library, leaf]
dependencies: []
reverse_dependencies: [anneal]
---

# LINK -- Optimization

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `link` |
| **Chapter** | 10 -- Optimization |
| **Purpose** | Linked-list utility for the simulated-annealing TSP solver. |
| **Status** | `library_only` |
| **Complexity** | `high` |
| **Fortran LOC** | 85 |
| **Subroutine** | `METROP` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/10_optimization/link/link.f` (85 lines)
- **Driver/demo:** _(none -- library-only or program-in-.f)_
- **Target:** `matarized/10_optimization/link/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    link
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `anneal` (10_optimization)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `ANS` | `LOGICAL` | (scalar) | parameter (input) | `double` |  |
| `DE` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `IORDER` | `INTEGER` | NCITY | local | `DFMatrixKokkos<int>(NCITY)` |  |
| `JDUM` | `INTEGER` | (scalar) | constant | `constexpr int JDUM = 1;` | constant = 1 |
| `JORDER` | `INTEGER` | MAXCIT | local | `DFMatrixKokkos<int>(MAXCIT)` |  |
| `MAXCIT` | `INTEGER` | (scalar) | constant | `constexpr int MAXCIT = 1000;` | constant = 1000 |
| `N` | `INTEGER` | 6 | local | `DFMatrixKokkos<int>(6)` |  |
| `T` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X` | `REAL` | NCITY | local | `DFMatrixKokkos<double>(NCITY)` |  |
| `XX` | `REAL` | 6 | local | `DFMatrixKokkos<double>(6)` |  |
| `Y` | `REAL` | NCITY | local | `DFMatrixKokkos<double>(NCITY)` |  |
| `YY` | `REAL` | 6 | local | `DFMatrixKokkos<double>(6)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,4

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K2: DO 11  J=1,NN

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K3: DO 11  J=1,6

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K4: DO 11  J=1,M1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K5: DO 12  J=1,M2

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K6: DO 13  J=1,M3

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.

### K7: DO 14  J=1,NCITY

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: 1, NN Array write(s) not indexed by loop variable: ALEN, IORDER, JORDER, N. Verify thread safety.


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
inline void link(double de, double t, double ans)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `METROP` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(link_matar_parallel CXX)

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


# --- Build the LINK example ---
add_executable(link main.cpp)
target_link_libraries(link matar_lib)

```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Reduction fusion:** If multiple reductions share the same loop bounds, consider fusing them into a single pass to reduce kernel launch overhead.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.
- **Hierarchical parallelism:** For deeply nested loops, consider `FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.

## 10. Validation Plan

### Reference Output

_No standalone Fortran driver (.dem) exists.  Validate by calling this routine from a dependent example._


### Serial Validation

```bash
cd matarized/10_optimization/link
mkdir -p build && cd build
cmake .. && make
./link > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/10_optimization/link/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/10_optimization/link
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./link > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./link > omp4_output.txt 2>&1
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
| **Conversion order** | 63 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | high (85 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `optimization`, `minimization`, `library`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
