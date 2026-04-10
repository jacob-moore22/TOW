---
example: bcuint
chapter: 03_interpolation
chapter_title: "Interpolation"
status: not_started
complexity: low
conversion_order: 127
priority: 0
tags: [interpolation, approximation]
dependencies: [bcucof]
reverse_dependencies: []
---

# BCUINT -- Interpolation

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `bcuint` |
| **Chapter** | 03 -- Interpolation |
| **Purpose** | Bicubic interpolation within a grid cell using bcucof. |
| **Status** | `not_started` |
| **Complexity** | `low` |
| **Fortran LOC** | 20 |
| **Subroutine** | `BCUINT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/03_interpolation/bcuint/bcuint.f` (20 lines)
- **Driver/demo:** `fortran/03_interpolation/bcuint/bcuint.dem`
- **Target:** `matarized/03_interpolation/bcuint/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `bcucof` (03_interpolation)

### Diagram

```mermaid
graph TD
    bcuint --> bcucof
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `ANS Y2` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ANSY` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ANSY1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `C` | `REAL` | 4, 4 | local | `DFMatrixKokkos<double>(4, 4)` |  |
| `X1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X1L` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X1U` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2L` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2U` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL` | 4 | parameter (input) | `DFMatrixKokkos<double>(4)` |  |
| `Y1` | `REAL` | 4 | parameter (input) | `DFMatrixKokkos<double>(4)` |  |
| `Y12` | `REAL` | 4 | parameter (input) | `DFMatrixKokkos<double>(4)` |  |
| `Y2` | `REAL` | 4 | parameter (input) | `DFMatrixKokkos<double>(4)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=4,1, step -1

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
inline void bcuint(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& y1, DFMatrixKokkos<double>& y2, DFMatrixKokkos<double>& y12, double x1l, double x1u, double x2l, double x2u, double x1, double x2, double ansy, double ansy1, double ans y2)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `BCUINT` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(bcuint_matar_parallel CXX)

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
set(INTERPOLATION_DIR    ${MATARIZED_ROOT}/03_interpolation)

# --- Build the BCUINT example ---
add_executable(bcuint main.cpp)
target_link_libraries(bcuint matar_lib)
target_include_directories(bcuint PRIVATE
    ${INTERPOLATION_DIR}/bcucof
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/03_interpolation/bcuint
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/03_interpolation/bcuint
mkdir -p build && cd build
cmake .. && make
./bcuint > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/03_interpolation/bcuint/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/03_interpolation/bcuint
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./bcuint > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./bcuint > omp4_output.txt 2>&1
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
| **Conversion order** | 127 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | low (20 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `bcucof` |
| **Tags** | `interpolation`, `approximation` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
