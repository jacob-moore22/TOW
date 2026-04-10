---
example: fixrts
chapter: 13_spectral_analysis
chapter_title: "Spectral Analysis"
status: not_started
complexity: medium
conversion_order: 91
priority: 1
tags: [spectral-analysis, signal-processing, cross-chapter]
dependencies: [laguer, zroots]
reverse_dependencies: [predic]
---

# FIXRTS -- Spectral Analysis

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `fixrts` |
| **Chapter** | 13 -- Spectral Analysis |
| **Purpose** | Stabilize a polynomial by moving roots inside the unit circle. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 31 |
| **Subroutine** | `FIXRTS` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/13_spectral_analysis/fixrts/fixrts.f` (31 lines)
- **Driver/demo:** `fortran/13_spectral_analysis/fixrts/fixrts.dem`
- **Target:** `matarized/13_spectral_analysis/fixrts/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `laguer` (09_root_finding)
  - `zroots` (09_root_finding)

### Diagram

```mermaid
graph TD
    fixrts --> laguer
    fixrts --> zroots
```

### Cross-Chapter Dependencies

- `laguer` from chapter 09
- `zroots` from chapter 09

## 4. Reverse Dependencies (examples that depend on this)

  - `predic` (13_spectral_analysis)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `COMPLEX` | NPMAX | local | `DFMatrixKokkos<double>(NPMAX)` |  |
| `D` | `REAL` | NPOLES | parameter (input) | `DFMatrixKokkos<double>(NPOLES)` |  |
| `NPMAX` | `INTEGER` | (scalar) | constant | `constexpr int NPMAX = 100;` | constant = 100 |
| `NPOLES` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `POLISH` | `LOGICAL` | (scalar) | local | `double` |  |
| `ROOTS` | `COMPLEX` | NPMAX | local | `DFMatrixKokkos<double>(NPMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=NPOLES,1, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K2: DO 12  J=1,NPOLES

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K3: DO 14  J=2,NPOLES

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K4: DO 13  I=J,2, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.

### K5: DO 15  J=1,NPOLES

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: A. Verify thread safety.


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
inline void fixrts(DFMatrixKokkos<double>& d, int npoles)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `FIXRTS` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(fixrts_matar_parallel CXX)

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
set(ROOTFINDING_DIR      ${MATARIZED_ROOT}/09_root_finding)

# --- Build the FIXRTS example ---
add_executable(fixrts main.cpp)
target_link_libraries(fixrts matar_lib)
target_include_directories(fixrts PRIVATE
    ${ROOTFINDING_DIR}/laguer
    ${ROOTFINDING_DIR}/zroots
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
cd fortran/13_spectral_analysis/fixrts
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/13_spectral_analysis/fixrts
mkdir -p build && cd build
cmake .. && make
./fixrts > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/13_spectral_analysis/fixrts/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/13_spectral_analysis/fixrts
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./fixrts > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./fixrts > omp4_output.txt 2>&1
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
| **Conversion order** | 91 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | medium (31 Fortran LOC, 2 dependencies) |
| **Prerequisite conversions** | `laguer`, `zroots` |
| **Tags** | `spectral-analysis`, `signal-processing`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
