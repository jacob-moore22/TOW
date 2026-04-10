---
example: cntab2
chapter: 14_statistics
chapter_title: "Statistics"
status: not_started
complexity: medium
conversion_order: 101
priority: 0
tags: [statistics, hypothesis-testing, leaf]
dependencies: []
reverse_dependencies: []
---

# CNTAB2 -- Statistics

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `cntab2` |
| **Chapter** | 14 -- Statistics |
| **Purpose** | Contingency table analysis using entropy-based measures. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 48 |
| **Subroutine** | `CNTAB2` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/14_statistics/cntab2/cntab2.f` (48 lines)
- **Driver/demo:** `fortran/14_statistics/cntab2/cntab2.dem`
- **Data files:** `TABLE.DAT`
- **Target:** `matarized/14_statistics/cntab2/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    cntab2
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `H` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HXGY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HYGX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `MAXI` | `INTEGER` | (scalar) | constant | `constexpr int MAXI = 100;` | constant = 100 |
| `MAXJ` | `INTEGER` | (scalar) | constant | `constexpr int MAXJ = 100;` | constant = 100 |
| `NI` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NJ` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NN` | `INTEGER` | NI, NJ | parameter (input) | `DFMatrixKokkos<int>(NI, NJ)` |  |
| `SUMI` | `REAL` | MAXI | local | `DFMatrixKokkos<double>(MAXI)` |  |
| `SUMJ` | `REAL` | MAXJ | local | `DFMatrixKokkos<double>(MAXJ)` |  |
| `TINY` | `REAL` | (scalar) | constant | `constexpr double TINY = 1.E-30;` | constant = 1.E-30 |
| `UXGY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `UXY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `UYGX` | `REAL` | (scalar) | parameter (input) | `double` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  I=1,NI

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K2: DO 11  J=1,NJ

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K3: DO 14  J=1,NJ

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K4: DO 13  I=1,NI

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K5: DO 15  I=1,NI

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K6: DO 16  J=1,NJ

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K7: DO 18  I=1,NI

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM

### K8: DO 17  J=1,NJ

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM


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
inline void cntab2(DFMatrixKokkos<double>& nn, int ni, int nj, double h, double hx, double hy, double hygx, double hxgy, double uygx, double uxgy, double uxy)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `CNTAB2` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(cntab2_matar_parallel CXX)

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


# --- Build the CNTAB2 example ---
add_executable(cntab2 main.cpp)
target_link_libraries(cntab2 matar_lib)

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
cd fortran/14_statistics/cntab2
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/14_statistics/cntab2
mkdir -p build && cd build
cmake .. && make
./cntab2 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/14_statistics/cntab2/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/14_statistics/cntab2
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./cntab2 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./cntab2 > omp4_output.txt 2>&1
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
| **Conversion order** | 101 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (48 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `statistics`, `hypothesis-testing`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
