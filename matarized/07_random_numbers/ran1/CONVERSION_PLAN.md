---
example: ran1
chapter: 07_random_numbers
chapter_title: "Random Numbers"
status: not_started
complexity: low
conversion_order: 3
priority: 25
tags: [random-number, stochastic, leaf]
dependencies: []
reverse_dependencies: [bnldev, expdev, gamdev, gasdev, poidev, mdian1, mdian2, smooft, avevar, chsone, chstwo, ftest, kendl1, ksone, kstwo, tptest, ttest, tutest, fit, lfit, medfit, mrqcof, mrqmin, rofunc, svdfit]
---

# RAN1 -- Random Numbers

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `ran1` |
| **Chapter** | 07 -- Random Numbers |
| **Purpose** | Minimal random number generator with Bays-Durham shuffle. |
| **Status** | `not_started` |
| **Complexity** | `low` |
| **Fortran LOC** | 30 |
| **Subroutine** | `RAN1` (function) |

## 2. Source Files

- **Fortran source:** `fortran/07_random_numbers/ran1/ran1.f` (30 lines)
- **Driver/demo:** `fortran/07_random_numbers/ran1/ran1.dem`
- **Target:** `matarized/07_random_numbers/ran1/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    ran1
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `bnldev` (07_random_numbers)
  - `expdev` (07_random_numbers)
  - `gamdev` (07_random_numbers)
  - `gasdev` (07_random_numbers)
  - `poidev` (07_random_numbers)
  - `mdian1` (08_sorting)
  - `mdian2` (08_sorting)
  - `smooft` (13_spectral_analysis)
  - `avevar` (14_statistics)
  - `chsone` (14_statistics)
  - `chstwo` (14_statistics)
  - `ftest` (14_statistics)
  - `kendl1` (14_statistics)
  - `ksone` (14_statistics)
  - `kstwo` (14_statistics)
  - `tptest` (14_statistics)
  - `ttest` (14_statistics)
  - `tutest` (14_statistics)
  - `fit` (15_curve_fitting)
  - `lfit` (15_curve_fitting)
  - `medfit` (15_curve_fitting)
  - `mrqcof` (15_curve_fitting)
  - `mrqmin` (15_curve_fitting)
  - `rofunc` (15_curve_fitting)
  - `svdfit` (15_curve_fitting)

> **Conversion note:** This routine is depended on by 25 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `IA1` | `INTEGER` | (scalar) | constant | `constexpr int IA1 = 7141;` | constant = 7141 |
| `IA2` | `INTEGER` | (scalar) | constant | `constexpr int IA2 = 8121;` | constant = 8121 |
| `IA3` | `INTEGER` | (scalar) | constant | `constexpr int IA3 = 4561;` | constant = 4561 |
| `IC1` | `INTEGER` | (scalar) | constant | `constexpr int IC1 = 54773;` | constant = 54773 |
| `IC2` | `INTEGER` | (scalar) | constant | `constexpr int IC2 = 28411;` | constant = 28411 |
| `IC3` | `INTEGER` | (scalar) | constant | `constexpr int IC3 = 51349;` | constant = 51349 |
| `IDUM` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `M1` | `INTEGER` | (scalar) | constant | `constexpr int M1 = 259200;` | constant = 259200 |
| `M2` | `INTEGER` | (scalar) | constant | `constexpr int M2 = 134456;` | constant = 134456 |
| `M3` | `INTEGER` | (scalar) | constant | `constexpr int M3 = 243000;` | constant = 243000 |
| `R` | `REAL` | 97 | local | `DFMatrixKokkos<double>(97)` |  |
| `RM1` | `REAL` | (scalar) | constant | `constexpr double RM1 = 3.8580247E-6;` | constant = 3.8580247E-6 |
| `RM2` | `REAL` | (scalar) | constant | `constexpr double RM2 = 7.4373773E-6;` | constant = 7.4373773E-6 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,97  (shuffle table initialization)

- **Thread safety:** `inherently_serial`
- **Recommended macro:** _serial `for` loop_
- **Notes:** The initialization loop fills a shuffle table using the LCG state. Each iteration depends on the LCG state from the previous iteration (IX1, IX2, IX3 are updated sequentially). **Parallelization note:** Random number generators are inherently sequential per-stream. For parallel use, create independent RNG streams per thread (e.g., different seeds or a counter-based RNG like Kokkos::Random). The Bays-Durham shuffle table makes per-call state dependencies even harder to break. **Strategy:** Use `Kokkos::Random_XorShift64_Pool` for parallel random number generation instead of porting ran1 directly.


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
inline double ran1(double idum)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `RAN1` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(ran1_matar_parallel CXX)

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


# --- Build the RAN1 example ---
add_executable(ran1 main.cpp)
target_link_libraries(ran1 matar_lib)

```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/07_random_numbers/ran1
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/07_random_numbers/ran1
mkdir -p build && cd build
cmake .. && make
./ran1 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/07_random_numbers/ran1/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/07_random_numbers/ran1
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./ran1 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./ran1 > omp4_output.txt 2>&1
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
| **Conversion order** | 3 of 202 |
| **Priority score** | 25 (reverse dependency count) |
| **Estimated effort** | low (30 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `random-number`, `stochastic`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
