---
example: eigsrt
chapter: 11_eigensystems
chapter_title: "Eigensystems"
status: not_started
complexity: low
conversion_order: 130
priority: 0
tags: [eigensystem, eigenvalue, matrix]
dependencies: [jacobi]
reverse_dependencies: []
---

# EIGSRT -- Eigensystems

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `eigsrt` |
| **Chapter** | 11 -- Eigensystems |
| **Purpose** | Sort eigenvalues and eigenvectors into descending order. |
| **Status** | `not_started` |
| **Complexity** | `low` |
| **Fortran LOC** | 24 |
| **Subroutine** | `EIGSRT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/11_eigensystems/eigsrt/eigsrt.f` (24 lines)
- **Driver/demo:** `fortran/11_eigensystems/eigsrt/eigsrt.dem`
- **Target:** `matarized/11_eigensystems/eigsrt/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `jacobi` (11_eigensystems)

### Diagram

```mermaid
graph TD
    eigsrt --> jacobi
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `D` | `REAL` | NP | parameter (input) | `DFMatrixKokkos<double>(NP)` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `V` | `REAL` | NP, NP | parameter (input) | `DFMatrixKokkos<double>(NP, NP)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 13  I=1,N-1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: D. Verify thread safety.

### K2: DO 11  J=I+1,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: D. Verify thread safety.

### K3: DO 12  J=1,N

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: D. Verify thread safety.


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
inline void eigsrt(DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& v, int n, int np)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `EIGSRT` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(eigsrt_matar_parallel CXX)

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
set(EIGENSYSTEMS_DIR     ${MATARIZED_ROOT}/11_eigensystems)

# --- Build the EIGSRT example ---
add_executable(eigsrt main.cpp)
target_link_libraries(eigsrt matar_lib)
target_include_directories(eigsrt PRIVATE
    ${EIGENSYSTEMS_DIR}/jacobi
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Loop ordering:** Verify innermost parallel index matches the fastest-varying array dimension for the chosen layout.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.
- **Hierarchical parallelism:** For deeply nested loops, consider `FOR_FIRST`/`FOR_SECOND` team-thread decomposition for better occupancy.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/11_eigensystems/eigsrt
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/11_eigensystems/eigsrt
mkdir -p build && cd build
cmake .. && make
./eigsrt > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/11_eigensystems/eigsrt/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/11_eigensystems/eigsrt
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./eigsrt > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./eigsrt > omp4_output.txt 2>&1
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
| **Conversion order** | 130 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | low (24 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `jacobi` |
| **Tags** | `eigensystem`, `eigenvalue`, `matrix` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
