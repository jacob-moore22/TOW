---
example: gaussj
chapter: 02_linear_algebra
chapter_title: "Linear Algebra"
status: not_started
complexity: high
conversion_order: 41
priority: 2
tags: [linear-algebra, matrix, leaf]
dependencies: []
reverse_dependencies: [lfit, mrqmin]
---

# GAUSSJ -- Linear Algebra

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `gaussj` |
| **Chapter** | 02 -- Linear Algebra |
| **Purpose** | Gauss-Jordan elimination for matrix inversion and solving linear systems Ax=b. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 72 |
| **Subroutine** | `GAUSSJ` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/02_linear_algebra/gaussj/gaussj.f` (72 lines)
- **Driver/demo:** `fortran/02_linear_algebra/gaussj/gaussj.dem`
- **Data files:** `MATRX1.DAT`
- **Target:** `matarized/02_linear_algebra/gaussj/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    gaussj
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `lfit` (15_curve_fitting)
  - `mrqmin` (15_curve_fitting)

> **Conversion note:** This routine is depended on by 2 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | NP, NP | parameter (input) | `DFMatrixKokkos<double>(NP, NP)` |  |
| `B` | `REAL` | NP, MP | parameter (input) | `DFMatrixKokkos<double>(NP, MP)` |  |
| `INDXC` | `INTEGER` | NMAX | local | `DFMatrixKokkos<int>(NMAX)` |  |
| `INDXR` | `INTEGER` | NMAX | local | `DFMatrixKokkos<int>(NMAX)` |  |
| `IPIV` | `INTEGER` | NMAX | local | `DFMatrixKokkos<int>(NMAX)` |  |
| `M` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `MP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 50;` | constant = 50 |
| `NP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 22  I=1,N  (OUTER ELIMINATION LOOP)

- **Thread safety:** `inherently_serial`
- **Recommended macro:** _serial `for` loop_
- **Notes:** This is the main Gauss-Jordan elimination loop. Each iteration selects a pivot (requiring knowledge of which columns are already pivoted), swaps rows, scales the pivot row, and eliminates entries in all other rows. Iterations have a hard sequential dependency (later pivots depend on the result of earlier eliminations). **Strategy:** Keep DO 22 as a serial `for` loop. Parallelize the *inner* loops (row swap DO 14/15, scaling DO 16/17, elimination DO 21 with nested DO 18/19) individually using `DO_ALL` with `MATAR_FENCE()` between dependent phases. The pivot search (DO 13 with nested DO 12) can use `DO_REDUCE_MAX` to find the largest element in parallel.

### K3: DO 13  J=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K4: DO 12  K=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K5: DO 14  L=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K6: DO 15  L=1,M

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K7: DO 16  L=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K8: DO 17  L=1,M

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K9: DO 21  LL=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K10: DO 18  L=1,N

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K11: DO 19  L=1,M

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K12: DO 24  L=N,1, step -1

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K13: DO 23  K=1,N

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
inline void gaussj(DFMatrixKokkos<double>& a, int n, int np, DFMatrixKokkos<double>& b, int m, int mp)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `GAUSSJ` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(gaussj_matar_parallel CXX)

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


# --- Build the GAUSSJ example ---
add_executable(gaussj main.cpp)
target_link_libraries(gaussj matar_lib)

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
cd fortran/02_linear_algebra/gaussj
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/02_linear_algebra/gaussj
mkdir -p build && cd build
cmake .. && make
./gaussj > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/02_linear_algebra/gaussj/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/02_linear_algebra/gaussj
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./gaussj > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./gaussj > omp4_output.txt 2>&1
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
| **Conversion order** | 41 of 202 |
| **Priority score** | 2 (reverse dependency count) |
| **Estimated effort** | high (72 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `linear-algebra`, `matrix`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
