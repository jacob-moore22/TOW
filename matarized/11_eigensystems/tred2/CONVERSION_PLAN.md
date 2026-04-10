---
example: tred2
chapter: 11_eigensystems
chapter_title: "Eigensystems"
status: not_started
complexity: high
conversion_order: 76
priority: 1
tags: [eigensystem, eigenvalue, matrix, leaf]
dependencies: []
reverse_dependencies: [tqli]
---

# TRED2 -- Eigensystems

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `tred2` |
| **Chapter** | 11 -- Eigensystems |
| **Purpose** | Householder reduction of a real symmetric matrix to tridiagonal form. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 81 |
| **Subroutine** | `TRED2` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/11_eigensystems/tred2/tred2.f` (81 lines)
- **Driver/demo:** `fortran/11_eigensystems/tred2/tred2.dem`
- **Target:** `matarized/11_eigensystems/tred2/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    tred2
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `tqli` (11_eigensystems)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | NP, NP | parameter (input) | `DFMatrixKokkos<double>(NP, NP)` |  |
| `D` | `REAL` | NP | parameter (input) | `DFMatrixKokkos<double>(NP)` |  |
| `E` | `REAL` | NP | parameter (input) | `DFMatrixKokkos<double>(NP)` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 18  I=N,2, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K2: DO 11  K=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K3: DO 12  K=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K4: DO 15  J=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K5: DO 13  K=1,J

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K6: DO 14  K=J+1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K7: DO 17  J=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K8: DO 16  K=1,J

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K9: DO 23  I=1,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K10: DO 21  J=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K11: DO 19  K=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K12: DO 20  K=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.

### K13: DO 22  J=1,L

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: F, G, H, SCALE Array write(s) not indexed by loop variable: D, E. Verify thread safety.


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
inline void tred2(DFMatrixKokkos<double>& a, int n, int np, DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& e)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `TRED2` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(tred2_matar_parallel CXX)

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


# --- Build the TRED2 example ---
add_executable(tred2 main.cpp)
target_link_libraries(tred2 matar_lib)

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
cd fortran/11_eigensystems/tred2
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/11_eigensystems/tred2
mkdir -p build && cd build
cmake .. && make
./tred2 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/11_eigensystems/tred2/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/11_eigensystems/tred2
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./tred2 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./tred2 > omp4_output.txt 2>&1
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
| **Conversion order** | 76 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | high (81 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `eigensystem`, `eigenvalue`, `matrix`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
