---
example: predic
chapter: 13_spectral_analysis
chapter_title: "Spectral Analysis"
status: not_started
complexity: high
conversion_order: 179
priority: 0
tags: [spectral-analysis, signal-processing, cross-chapter]
dependencies: [laguer, zroots, fixrts, memcof]
reverse_dependencies: []
---

# PREDIC -- Spectral Analysis

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `predic` |
| **Chapter** | 13 -- Spectral Analysis |
| **Purpose** | Linear prediction of future data using maximum entropy coefficients. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 21 |
| **Subroutine** | `PREDIC` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/13_spectral_analysis/predic/predic.f` (21 lines)
- **Driver/demo:** `fortran/13_spectral_analysis/predic/predic.dem`
- **Target:** `matarized/13_spectral_analysis/predic/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `laguer` (09_root_finding)
  - `zroots` (09_root_finding)
  - `fixrts` (13_spectral_analysis)
  - `memcof` (13_spectral_analysis)

### Diagram

```mermaid
graph TD
    predic --> laguer
    predic --> zroots
    predic --> fixrts
    predic --> memcof
```

### Cross-Chapter Dependencies

- `laguer` from chapter 09
- `zroots` from chapter 09

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `D` | `REAL` | NPOLES | parameter (input) | `DFMatrixKokkos<double>(NPOLES)` |  |
| `DATA` | `REAL` | NDATA | parameter (input) | `DFMatrixKokkos<double>(NDATA)` |  |
| `FUTURE` | `REAL` | NFUT | parameter (input) | `DFMatrixKokkos<double>(NFUT)` |  |
| `NDATA` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NFUT` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NPMAX` | `INTEGER` | (scalar) | constant | `constexpr int NPMAX = 100;` | constant = 100 |
| `NPOLES` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `REG` | `REAL` | NPMAX | local | `DFMatrixKokkos<double>(NPMAX)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,NPOLES

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM Array write(s) not indexed by loop variable: REG. Verify thread safety.

### K2: DO 14  J=1,NFUT

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM Array write(s) not indexed by loop variable: REG. Verify thread safety.

### K3: DO 12  K=1,NPOLES

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM Array write(s) not indexed by loop variable: REG. Verify thread safety.

### K4: DO 13  K=NPOLES,2, step -1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: SUM Array write(s) not indexed by loop variable: REG. Verify thread safety.


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
inline void predic(DFMatrixKokkos<double>& data, int ndata, DFMatrixKokkos<double>& d, int npoles, DFMatrixKokkos<double>& future, int nfut)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `PREDIC` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(predic_matar_parallel CXX)

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
set(SPECTRALANALYSIS_DIR ${MATARIZED_ROOT}/13_spectral_analysis)

# --- Build the PREDIC example ---
add_executable(predic main.cpp)
target_link_libraries(predic matar_lib)
target_include_directories(predic PRIVATE
    ${ROOTFINDING_DIR}/laguer
    ${ROOTFINDING_DIR}/zroots
    ${SPECTRALANALYSIS_DIR}/fixrts
    ${SPECTRALANALYSIS_DIR}/memcof
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
cd fortran/13_spectral_analysis/predic
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/13_spectral_analysis/predic
mkdir -p build && cd build
cmake .. && make
./predic > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/13_spectral_analysis/predic/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/13_spectral_analysis/predic
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./predic > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./predic > omp4_output.txt 2>&1
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
| **Conversion order** | 179 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (21 Fortran LOC, 4 dependencies) |
| **Prerequisite conversions** | `laguer`, `zroots`, `fixrts`, `memcof` |
| **Tags** | `spectral-analysis`, `signal-processing`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
