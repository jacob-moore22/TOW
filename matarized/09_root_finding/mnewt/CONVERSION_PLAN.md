---
example: mnewt
chapter: 09_root_finding
chapter_title: "Root Finding"
status: not_started
complexity: high
conversion_order: 178
priority: 0
tags: [root-finding, nonlinear-equation, cross-chapter]
dependencies: [lubksb, ludcmp, beta, gammln]
reverse_dependencies: []
---

# MNEWT -- Root Finding

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `mnewt` |
| **Chapter** | 09 -- Root Finding |
| **Purpose** | Multi-dimensional Newton-Raphson root finding for nonlinear systems. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 22 |
| **Subroutine** | `MNEWT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/09_root_finding/mnewt/mnewt.f` (22 lines)
- **Driver/demo:** `fortran/09_root_finding/mnewt/mnewt.dem`
- **Target:** `matarized/09_root_finding/mnewt/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `lubksb` (02_linear_algebra)
  - `ludcmp` (02_linear_algebra)
  - `beta` (06_special_functions)
  - `gammln` (06_special_functions)

### Diagram

```mermaid
graph TD
    mnewt --> lubksb
    mnewt --> ludcmp
    mnewt --> beta
    mnewt --> gammln
```

### Cross-Chapter Dependencies

- `lubksb` from chapter 02
- `ludcmp` from chapter 02
- `beta` from chapter 06
- `gammln` from chapter 06

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `ALPHA` | `REAL` | NP, NP | local | `DFMatrixKokkos<double>(NP, NP)` |  |
| `BETA` | `REAL` | NP | local | `DFMatrixKokkos<double>(NP)` |  |
| `INDX` | `INTEGER` | NP | local | `DFMatrixKokkos<int>(NP)` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NP` | `INTEGER` | (scalar) | constant | `constexpr int NP = 15;` | constant = 15 |
| `NTRIAL` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `TOLF` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `TOLX` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X` | `REAL` | NP | parameter (input) | `DFMatrixKokkos<double>(NP)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 13  K=1,NTRIAL

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: ERRF, ERRX

### K2: DO 11  I=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: ERRF, ERRX

### K3: DO 12  I=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: ERRF, ERRX


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
inline void mnewt(int ntrial, DFMatrixKokkos<double>& x, int n, double tolx, double tolf)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `MNEWT` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(mnewt_matar_parallel CXX)

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
set(LINEARALGEBRA_DIR    ${MATARIZED_ROOT}/02_linear_algebra)
set(SPECIALFUNCTIONS_DIR ${MATARIZED_ROOT}/06_special_functions)

# --- Build the MNEWT example ---
add_executable(mnewt main.cpp)
target_link_libraries(mnewt matar_lib)
target_include_directories(mnewt PRIVATE
    ${LINEARALGEBRA_DIR}/lubksb
    ${LINEARALGEBRA_DIR}/ludcmp
    ${SPECIALFUNCTIONS_DIR}/beta
    ${SPECIALFUNCTIONS_DIR}/gammln
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
cd fortran/09_root_finding/mnewt
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/09_root_finding/mnewt
mkdir -p build && cd build
cmake .. && make
./mnewt > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/09_root_finding/mnewt/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/09_root_finding/mnewt
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./mnewt > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./mnewt > omp4_output.txt 2>&1
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
| **Conversion order** | 178 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (22 Fortran LOC, 4 dependencies) |
| **Prerequisite conversions** | `lubksb`, `ludcmp`, `beta`, `gammln` |
| **Tags** | `root-finding`, `nonlinear-equation`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
