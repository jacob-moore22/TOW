---
example: splin2
chapter: 03_interpolation
chapter_title: "Interpolation"
status: not_started
complexity: medium
conversion_order: 172
priority: 0
tags: [interpolation, approximation]
dependencies: [splie2, spline, splint]
reverse_dependencies: []
---

# SPLIN2 -- Interpolation

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `splin2` |
| **Chapter** | 03 -- Interpolation |
| **Purpose** | Evaluate a 2D bicubic spline interpolation at a given point. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 16 |
| **Subroutine** | `SPLIN2` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/03_interpolation/splin2/splin2.f` (16 lines)
- **Driver/demo:** `fortran/03_interpolation/splin2/splin2.dem`
- **Target:** `matarized/03_interpolation/splin2/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `splie2` (03_interpolation)
  - `spline` (03_interpolation)
  - `splint` (03_interpolation)

### Diagram

```mermaid
graph TD
    splin2 --> splie2
    splin2 --> spline
    splin2 --> splint
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `M` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NN` | `INTEGER` | (scalar) | constant | `constexpr int NN = 100;` | constant = 100 |
| `X1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X1A` | `REAL` | M | parameter (input) | `DFMatrixKokkos<double>(M)` |  |
| `X2` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2A` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `Y` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y2A` | `REAL` | M, N | parameter (input) | `DFMatrixKokkos<double>(M, N)` |  |
| `Y2TMP` | `REAL` | NN | local | `DFMatrixKokkos<double>(NN)` |  |
| `YA` | `REAL` | M, N | parameter (input) | `DFMatrixKokkos<double>(M, N)` |  |
| `YTMP` | `REAL` | NN | local | `DFMatrixKokkos<double>(NN)` |  |
| `YYTMP` | `REAL` | NN | local | `DFMatrixKokkos<double>(NN)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  J=1,M

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 11  K=1,N

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
inline void splin2(DFMatrixKokkos<double>& x1a, DFMatrixKokkos<double>& x2a, DFMatrixKokkos<double>& ya, DFMatrixKokkos<double>& y2a, int m, int n, double x1, double x2, double y)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `SPLIN2` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(splin2_matar_parallel CXX)

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

# --- Build the SPLIN2 example ---
add_executable(splin2 main.cpp)
target_link_libraries(splin2 matar_lib)
target_include_directories(splin2 PRIVATE
    ${INTERPOLATION_DIR}/splie2
    ${INTERPOLATION_DIR}/spline
    ${INTERPOLATION_DIR}/splint
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
cd fortran/03_interpolation/splin2
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/03_interpolation/splin2
mkdir -p build && cd build
cmake .. && make
./splin2 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/03_interpolation/splin2/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/03_interpolation/splin2
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./splin2 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./splin2 > omp4_output.txt 2>&1
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
| **Conversion order** | 172 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (16 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `splie2`, `spline`, `splint` |
| **Tags** | `interpolation`, `approximation` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
