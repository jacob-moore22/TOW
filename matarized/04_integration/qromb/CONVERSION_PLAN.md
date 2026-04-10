---
example: qromb
chapter: 04_integration
chapter_title: "Integration"
status: not_started
complexity: medium
conversion_order: 152
priority: 0
tags: [integration, quadrature, cross-chapter]
dependencies: [polint, trapzd]
reverse_dependencies: []
---

# QROMB -- Integration

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `qromb` |
| **Chapter** | 04 -- Integration |
| **Purpose** | Romberg integration using Richardson extrapolation on the trapezoidal rule. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 17 |
| **Subroutine** | `QROMB` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/04_integration/qromb/qromb.f` (17 lines)
- **Driver/demo:** `fortran/04_integration/qromb/qromb.dem`
- **Target:** `matarized/04_integration/qromb/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `polint` (03_interpolation)
  - `trapzd` (04_integration)

### Diagram

```mermaid
graph TD
    qromb --> polint
    qromb --> trapzd
```

### Cross-Chapter Dependencies

- `polint` from chapter 03

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `B` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `EPS` | `REAL` | (scalar) | constant | `constexpr double EPS = 1.E-6;` | constant = 1.E-6 |
| `FUNC` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `H` | `REAL` | JMAXP | local | `DFMatrixKokkos<double>(JMAXP)` |  |
| `JMAX` | `INTEGER` | (scalar) | constant | `constexpr int JMAX = 20;` | constant = 20 |
| `JMAXP` | `INTEGER` | (scalar) | constant | `constexpr int JMAXP = JMAX+1;` | constant = JMAX+1 |
| `K` | `INTEGER` | (scalar) | constant | `constexpr int K = 5;` | constant = 5 |
| `KM` | `INTEGER` | (scalar) | constant | `constexpr int KM = 4;` | constant = 4 |
| `S` | `REAL` | JMAXP | local | `DFMatrixKokkos<double>(JMAXP)` |  |
| `SS` | `REAL` | (scalar) | parameter (input) | `double` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,JMAX

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: H. Verify thread safety.


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
inline void qromb(double func, double a, double b, double ss)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `QROMB` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(qromb_matar_parallel CXX)

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
set(INTEGRATION_DIR      ${MATARIZED_ROOT}/04_integration)

# --- Build the QROMB example ---
add_executable(qromb main.cpp)
target_link_libraries(qromb matar_lib)
target_include_directories(qromb PRIVATE
    ${INTERPOLATION_DIR}/polint
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
cd fortran/04_integration/qromb
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/04_integration/qromb
mkdir -p build && cd build
cmake .. && make
./qromb > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/04_integration/qromb/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/04_integration/qromb
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./qromb > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./qromb > omp4_output.txt 2>&1
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
| **Conversion order** | 152 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (17 Fortran LOC, 2 dependencies) |
| **Prerequisite conversions** | `polint`, `trapzd` |
| **Tags** | `integration`, `quadrature`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
