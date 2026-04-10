---
example: des
chapter: 07_random_numbers
chapter_title: "Random Numbers"
status: not_started
complexity: medium
conversion_order: 50
priority: 2
tags: [random-number, stochastic]
dependencies: [desks]
reverse_dependencies: [ran4, kendl1]
---

# DES -- Random Numbers

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `des` |
| **Chapter** | 07 -- Random Numbers |
| **Purpose** | Data Encryption Standard (DES) block cipher implementation. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 41 |
| **Subroutine** | `DES` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/07_random_numbers/des/des.f` (41 lines)
- **Driver/demo:** `fortran/07_random_numbers/des/des.dem`
- **Data files:** `DESTST.DAT`
- **Target:** `matarized/07_random_numbers/des/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `desks` (07_random_numbers)

### Diagram

```mermaid
graph TD
    des --> desks
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `ran4` (07_random_numbers)
  - `kendl1` (14_statistics)

> **Conversion note:** This routine is depended on by 2 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `ICF` | `INTEGER` | 32 | local | `DFMatrixKokkos<int>(32)` |  |
| `INPUT` | `INTEGER` | 64 | parameter (input) | `DFMatrixKokkos<int>(64)` |  |
| `IP` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `IPM` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `ISW` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `ITMP` | `INTEGER` | 64 | local | `DFMatrixKokkos<int>(64)` |  |
| `JOTPUT` | `INTEGER` | 64 | parameter (input) | `DFMatrixKokkos<int>(64)` |  |
| `KEY` | `INTEGER` | 64 | parameter (input) | `DFMatrixKokkos<int>(64)` |  |
| `KNS` | `INTEGER` | 48, 16 | local | `DFMatrixKokkos<int>(48, 16)` |  |
| `NEWKEY` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=1,16

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K2: DO 12  J=1,64

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K3: DO 14  I=1,16

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** None

### K4: DO 13  J=1,32

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K5: DO 15  J=1,32

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.

### K6: DO 16  J=1,64

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Stencil access pattern detected -- verify neighbor independence.


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
inline void des(DFMatrixKokkos<double>& input, DFMatrixKokkos<double>& key, int newkey, int isw, DFMatrixKokkos<double>& jotput)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `DES` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(des_matar_parallel CXX)

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
set(RANDOMNUMBERS_DIR    ${MATARIZED_ROOT}/07_random_numbers)

# --- Build the DES example ---
add_executable(des main.cpp)
target_link_libraries(des matar_lib)
target_include_directories(des PRIVATE
    ${RANDOMNUMBERS_DIR}/desks
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
cd fortran/07_random_numbers/des
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/07_random_numbers/des
mkdir -p build && cd build
cmake .. && make
./des > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/07_random_numbers/des/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/07_random_numbers/des
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./des > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./des > omp4_output.txt 2>&1
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
| **Conversion order** | 50 of 202 |
| **Priority score** | 2 (reverse dependency count) |
| **Estimated effort** | medium (41 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `desks` |
| **Tags** | `random-number`, `stochastic` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
