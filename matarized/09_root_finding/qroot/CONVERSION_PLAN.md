---
example: qroot
chapter: 09_root_finding
chapter_title: "Root Finding"
status: not_started
complexity: medium
conversion_order: 135
priority: 0
tags: [root-finding, nonlinear-equation, cross-chapter]
dependencies: [poldiv]
reverse_dependencies: []
---

# QROOT -- Root Finding

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `qroot` |
| **Chapter** | 09 -- Root Finding |
| **Purpose** | Find a quadratic factor of a polynomial by Bairstow's method. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 32 |
| **Subroutine** | `QROOT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/09_root_finding/qroot/qroot.f` (32 lines)
- **Driver/demo:** `fortran/09_root_finding/qroot/qroot.dem`
- **Target:** `matarized/09_root_finding/qroot/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `poldiv` (05_evaluation_of_functions)

### Diagram

```mermaid
graph TD
    qroot --> poldiv
```

### Cross-Chapter Dependencies

- `poldiv` from chapter 05

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `B` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `C` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `D` | `REAL` | 3 | local | `DFMatrixKokkos<double>(3)` |  |
| `EPS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `ITMAX` | `INTEGER` | (scalar) | constant | `constexpr int ITMAX = 20;` | constant = 20 |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 20;` | constant = 20 |
| `P` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `Q` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `QQ` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `REM` | `REAL` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `TINY` | `REAL` | (scalar) | constant | `constexpr double TINY = 1.0E-6;` | constant = 1.0E-6 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  ITER=1,ITMAX

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: B, C Array write(s) not indexed by loop variable: D, Q. Verify thread safety.

### K2: DO 11  I=N-1,1, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: B, C Array write(s) not indexed by loop variable: D, Q. Verify thread safety.


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
inline void qroot(DFMatrixKokkos<double>& p, int n, double b, double c, double eps)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `QROOT` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(qroot_matar_parallel CXX)

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
set(EVALUATIONOFFUNCTIONS_DIR ${MATARIZED_ROOT}/05_evaluation_of_functions)

# --- Build the QROOT example ---
add_executable(qroot main.cpp)
target_link_libraries(qroot matar_lib)
target_include_directories(qroot PRIVATE
    ${EVALUATIONOFFUNCTIONS_DIR}/poldiv
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Reduction fusion:** If multiple reductions share the same loop bounds, consider fusing them into a single pass to reduce kernel launch overhead.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/09_root_finding/qroot
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/09_root_finding/qroot
mkdir -p build && cd build
cmake .. && make
./qroot > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/09_root_finding/qroot/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/09_root_finding/qroot
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./qroot > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./qroot > omp4_output.txt 2>&1
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
| **Conversion order** | 135 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (32 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `poldiv` |
| **Tags** | `root-finding`, `nonlinear-equation`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
