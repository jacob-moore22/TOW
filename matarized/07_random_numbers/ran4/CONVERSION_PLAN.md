---
example: ran4
chapter: 07_random_numbers
chapter_title: "Random Numbers"
status: not_started
complexity: medium
conversion_order: 95
priority: 1
tags: [random-number, stochastic, cross-chapter]
dependencies: [des, desks, avevar]
reverse_dependencies: [kendl1]
---

# RAN4 -- Random Numbers

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `ran4` |
| **Chapter** | 07 -- Random Numbers |
| **Purpose** | Random number generator using the Data Encryption Standard. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 34 |
| **Subroutine** | `RAN4` (function) |

## 2. Source Files

- **Fortran source:** `fortran/07_random_numbers/ran4/ran4.f` (34 lines)
- **Driver/demo:** `fortran/07_random_numbers/ran4/ran4.dem`
- **Target:** `matarized/07_random_numbers/ran4/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `des` (07_random_numbers)
  - `desks` (07_random_numbers)
  - `avevar` (14_statistics)

### Diagram

```mermaid
graph TD
    ran4 --> des
    ran4 --> desks
    ran4 --> avevar
```

### Cross-Chapter Dependencies

- `avevar` from chapter 14

## 4. Reverse Dependencies (examples that depend on this)

  - `kendl1` (14_statistics)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `IA` | `INTEGER` | (scalar) | constant | `constexpr int IA = 430;` | constant = 430 |
| `IC` | `INTEGER` | (scalar) | constant | `constexpr int IC = 2531;` | constant = 2531 |
| `IDUM` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `IM` | `INTEGER` | (scalar) | constant | `constexpr int IM = 11979;` | constant = 11979 |
| `INP` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `JOT` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `KEY` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `NACC` | `INTEGER` | (scalar) | constant | `constexpr int NACC = 24;` | constant = 24 |
| `POW` | `REAL` | 65 | local | `DFMatrixKokkos<double>(65)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,64

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: RAN4 Array write(s) not indexed by loop variable: INP, POW. Verify thread safety.

### K2: DO 12  J=64,2, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: RAN4 Array write(s) not indexed by loop variable: INP, POW. Verify thread safety.

### K3: DO 13  J=1,NACC

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: RAN4 Array write(s) not indexed by loop variable: INP, POW. Verify thread safety.


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
inline double ran4(double idum)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `RAN4` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(ran4_matar_parallel CXX)

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
set(RANDOMNUMBERS_DIR    ${MATARIZED_ROOT}/07_random_numbers)
set(STATISTICS_DIR       ${MATARIZED_ROOT}/14_statistics)

# --- Build the RAN4 example ---
add_executable(ran4 main.cpp)
target_link_libraries(ran4 matar_lib)
target_include_directories(ran4 PRIVATE
    ${RANDOMNUMBERS_DIR}/des
    ${RANDOMNUMBERS_DIR}/desks
    ${STATISTICS_DIR}/avevar
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
cd fortran/07_random_numbers/ran4
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/07_random_numbers/ran4
mkdir -p build && cd build
cmake .. && make
./ran4 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/07_random_numbers/ran4/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/07_random_numbers/ran4
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./ran4 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./ran4 > omp4_output.txt 2>&1
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
| **Conversion order** | 95 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | medium (34 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `des`, `desks`, `avevar` |
| **Tags** | `random-number`, `stochastic`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
