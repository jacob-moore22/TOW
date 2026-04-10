---
example: bessi0
chapter: 06_special_functions
chapter_title: "Special Functions"
status: not_started
complexity: low
conversion_order: 28
priority: 3
tags: [special-function, mathematical-function, leaf]
dependencies: []
reverse_dependencies: [bessi, bessk, bessk0]
---

# BESSI0 -- Special Functions

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `bessi0` |
| **Chapter** | 06 -- Special Functions |
| **Purpose** | Modified Bessel function I0(x) using polynomial approximation. |
| **Status** | `not_started` |
| **Complexity** | `low` |
| **Fortran LOC** | 21 |
| **Subroutine** | `BESSI0` (function) |

## 2. Source Files

- **Fortran source:** `fortran/06_special_functions/bessi0/bessi0.f` (21 lines)
- **Driver/demo:** `fortran/06_special_functions/bessi0/bessi0.dem`
- **Data files:** `FNCVAL.DAT`
- **Target:** `matarized/06_special_functions/bessi0/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  (none)

### Diagram

```mermaid
graph TD
    bessi0
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `bessi` (06_special_functions)
  - `bessk` (06_special_functions)
  - `bessk0` (06_special_functions)

> **Conversion note:** This routine is depended on by 3 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `P1` | `REAL*8` | (scalar) | local | `double` |  |
| `P2` | `REAL*8` | (scalar) | local | `double` |  |
| `P3` | `REAL*8` | (scalar) | local | `double` |  |
| `P4` | `REAL*8` | (scalar) | local | `double` |  |
| `P5` | `REAL*8` | (scalar) | local | `double` |  |
| `P6` | `REAL*8` | (scalar) | local | `double` |  |
| `P7` | `REAL*8` | (scalar) | local | `double` |  |
| `Q1` | `REAL*8` | (scalar) | local | `double` |  |
| `Q2` | `REAL*8` | (scalar) | local | `double` |  |
| `Q3` | `REAL*8` | (scalar) | local | `double` |  |
| `Q4` | `REAL*8` | (scalar) | local | `double` |  |
| `Q5` | `REAL*8` | (scalar) | local | `double` |  |
| `Q6` | `REAL*8` | (scalar) | local | `double` |  |
| `Q7` | `REAL*8` | (scalar) | local | `double` |  |
| `Q8` | `REAL*8` | (scalar) | local | `double` |  |
| `Q9` | `REAL*8` | (scalar) | local | `double` |  |
| `X` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL*8` | (scalar) | local | `double` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

_No DO loops detected in source._

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
inline double bessi0(double x)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `BESSI0` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(bessi0_matar_parallel CXX)

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


# --- Build the BESSI0 example ---
add_executable(bessi0 main.cpp)
target_link_libraries(bessi0 matar_lib)

```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/06_special_functions/bessi0
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/06_special_functions/bessi0
mkdir -p build && cd build
cmake .. && make
./bessi0 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/06_special_functions/bessi0/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/06_special_functions/bessi0
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./bessi0 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./bessi0 > omp4_output.txt 2>&1
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
| **Conversion order** | 28 of 202 |
| **Priority score** | 3 (reverse dependency count) |
| **Estimated effort** | low (21 Fortran LOC, 0 dependencies) |
| **Prerequisite conversions** | (none -- leaf node) |
| **Tags** | `special-function`, `mathematical-function`, `leaf` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
