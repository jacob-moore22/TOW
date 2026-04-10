---
example: ludcmp
chapter: 02_linear_algebra
chapter_title: "Linear Algebra"
status: not_started
complexity: high
conversion_order: 17
priority: 5
tags: [linear-algebra, matrix, leaf]
dependencies: []
reverse_dependencies: [lubksb, mprove, mnewt, shoot, shootf]
---

# LUDCMP -- Linear Algebra

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `ludcmp` |
| **Chapter** | 02 -- Linear Algebra |
| **Purpose** | LU decomposition of a square matrix using Crout's algorithm with partial pivoting. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 61 |
| **Subroutine** | `LUDCMP` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/02_linear_algebra/ludcmp/ludcmp.f` (61 lines)
- **Driver/demo:** `fortran/02_linear_algebra/ludcmp/ludcmp.dem`
- **Data files:** `MATRX1.DAT`
- **Target:** `matarized/02_linear_algebra/ludcmp/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    ludcmp
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `lubksb` (02_linear_algebra)
  - `mprove` (02_linear_algebra)
  - `mnewt` (09_root_finding)
  - `shoot` (17_boundary_value_problems)
  - `shootf` (17_boundary_value_problems)

> **Conversion note:** This routine is depended on by 5 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | NP, NP | parameter (input) | `DFMatrixKokkos<double>(NP, NP)` |  |
| `D` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `INDX` | `INTEGER` | N | parameter (input) | `DFMatrixKokkos<int>(N)` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 100;` | constant = 100 |
| `NP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `TINY` | `REAL` | (scalar) | constant | `constexpr double TINY = 1.0E-20;` | constant = 1.0E-20 |
| `VV` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12 I / DO 11 J  (row scaling: find max abs per row)

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL` for outer loop, `DO_REDUCE_MAX` for inner loop
- **Notes:** The outer loop iterates over rows; the inner loop finds the max absolute value per row. Rows are independent. Each row's max can be found with a parallel reduction.

### K2: DO 19  J=1,N  (OUTER CROUT DECOMPOSITION LOOP)

- **Thread safety:** `inherently_serial`
- **Recommended macro:** _serial `for` loop_
- **Notes:** This is the column-by-column Crout decomposition. Each column J depends on the results of all prior columns (triangular solve dependency). Must remain serial. Contains inner loops K3-K5 that can be parallelized individually.

### K3: DO 14 I / DO 13 K  (upper triangle: compute U elements)

- **Thread safety:** The outer DO 14 loop is parallelizable (each row I independent within the same column J). Inner DO 13 is a reduction for each row.
- **Recommended macro:** `DO_ALL` for I loop, inner K loop as serial accumulation per thread
- **Notes:** `SUM = A(I,K)*A(K,J)` accumulated over K is a dot product per row. Each row can be computed in parallel. For small matrices, keep inner K serial per thread.

### K4: DO 16 I / DO 15 K  (lower triangle: compute L elements + find pivot)

- **Thread safety:** `DO_ALL` for I loop + `DO_REDUCE_MAX` for pivot search
- **Recommended macro:** `DO_ALL` with per-thread dot product, then `DO_REDUCE_MAX` for pivot
- **Notes:** Similar structure to K3 but for the lower triangle. The pivot search (AAMAX) across rows is a max-reduction. The row swap after pivoting is a separate parallel kernel.

### K6: DO 16  I=J,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K7: DO 15  K=1,J-1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K8: DO 17  K=1,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K9: DO 18  I=J+1,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.


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
inline void ludcmp(DFMatrixKokkos<double>& a, int n, int np, DFMatrixKokkos<double>& indx, double d)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `LUDCMP` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(ludcmp_matar_parallel CXX)

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


# --- Build the LUDCMP example ---
add_executable(ludcmp main.cpp)
target_link_libraries(ludcmp matar_lib)

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
cd fortran/02_linear_algebra/ludcmp
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/02_linear_algebra/ludcmp
mkdir -p build && cd build
cmake .. && make
./ludcmp > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/02_linear_algebra/ludcmp/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/02_linear_algebra/ludcmp
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./ludcmp > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./ludcmp > omp4_output.txt 2>&1
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
| **Conversion order** | 17 of 202 |
| **Priority score** | 5 (reverse dependency count) |
| **Estimated effort** | high (61 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `linear-algebra`, `matrix`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
