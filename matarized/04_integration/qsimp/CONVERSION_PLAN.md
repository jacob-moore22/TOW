---
example: qsimp
chapter: 04_integration
chapter_title: "Integration"
status: not_started
complexity: low
conversion_order: 52
priority: 2
tags: [integration, quadrature]
dependencies: [trapzd]
reverse_dependencies: [cel, el2]
---

# QSIMP -- Integration

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `qsimp` |
| **Chapter** | 04 -- Integration |
| **Purpose** | Integration using Simpson's rule to convergence. |
| **Status** | `not_started` |
| **Complexity** | `low` |
| **Fortran LOC** | 14 |
| **Subroutine** | `QSIMP` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/04_integration/qsimp/qsimp.f` (14 lines)
- **Driver/demo:** `fortran/04_integration/qsimp/qsimp.dem`
- **Target:** `matarized/04_integration/qsimp/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `trapzd` (04_integration)

### Diagram

```mermaid
graph TD
    qsimp --> trapzd
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `cel` (06_special_functions)
  - `el2` (06_special_functions)

> **Conversion note:** This routine is depended on by 2 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `B` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `EPS` | `REAL` | (scalar) | constant | `constexpr double EPS = 1.E-6;` | constant = 1.E-6 |
| `FUNC` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `JMAX` | `INTEGER` | (scalar) | constant | `constexpr int JMAX = 20;` | constant = 20 |
| `S` | `REAL` | (scalar) | parameter (input) | `double` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,JMAX

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
inline void qsimp(double func, double a, double b, double s)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `QSIMP` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(qsimp_matar_parallel CXX)

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
set(INTEGRATION_DIR      ${MATARIZED_ROOT}/04_integration)

# --- Build the QSIMP example ---
add_executable(qsimp main.cpp)
target_link_libraries(qsimp matar_lib)
target_include_directories(qsimp PRIVATE
    ${INTEGRATION_DIR}/trapzd
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
cd fortran/04_integration/qsimp
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/04_integration/qsimp
mkdir -p build && cd build
cmake .. && make
./qsimp > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/04_integration/qsimp/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/04_integration/qsimp
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./qsimp > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./qsimp > omp4_output.txt 2>&1
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
| **Conversion order** | 52 of 202 |
| **Priority score** | 2 (reverse dependency count) |
| **Estimated effort** | low (14 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `trapzd` |
| **Tags** | `integration`, `quadrature` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
