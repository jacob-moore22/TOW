---
example: correl
chapter: 13_spectral_analysis
chapter_title: "Spectral Analysis"
status: not_started
complexity: medium
conversion_order: 162
priority: 0
tags: [spectral-analysis, signal-processing, cross-chapter]
dependencies: [four1, realft, twofft]
reverse_dependencies: []
---

# CORREL -- Spectral Analysis

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `correl` |
| **Chapter** | 13 -- Spectral Analysis |
| **Purpose** | FFT-based cross-correlation of two real data sets. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 14 |
| **Subroutine** | `CORREL` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/13_spectral_analysis/correl/correl.f` (14 lines)
- **Driver/demo:** `fortran/13_spectral_analysis/correl/correl.dem`
- **Target:** `matarized/13_spectral_analysis/correl/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `four1` (12_fourier_transform)
  - `realft` (12_fourier_transform)
  - `twofft` (12_fourier_transform)

### Diagram

```mermaid
graph TD
    correl --> four1
    correl --> realft
    correl --> twofft
```

### Cross-Chapter Dependencies

- `four1` from chapter 12
- `realft` from chapter 12
- `twofft` from chapter 12

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `ANS` | `COMPLEX` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `DATA1` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `DATA2` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `FFT` | `COMPLEX` | NMAX | local | `DFMatrixKokkos<double>(NMAX)` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NMAX` | `INTEGER` | (scalar) | constant | `constexpr int NMAX = 8192;` | constant = 8192 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  I=1,N/2+1

- **Thread safety:** `unsafe_review`
- **Recommended macro:** `DO_ALL`
- **Notes:** Array write(s) not indexed by loop variable: ANS. Verify thread safety.


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
inline void correl(DFMatrixKokkos<double>& data1, DFMatrixKokkos<double>& data2, int n, DFMatrixKokkos<double>& ans)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `CORREL` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(correl_matar_parallel CXX)

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
set(FOURIERTRANSFORM_DIR ${MATARIZED_ROOT}/12_fourier_transform)

# --- Build the CORREL example ---
add_executable(correl main.cpp)
target_link_libraries(correl matar_lib)
target_include_directories(correl PRIVATE
    ${FOURIERTRANSFORM_DIR}/four1
    ${FOURIERTRANSFORM_DIR}/realft
    ${FOURIERTRANSFORM_DIR}/twofft
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
cd fortran/13_spectral_analysis/correl
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/13_spectral_analysis/correl
mkdir -p build && cd build
cmake .. && make
./correl > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/13_spectral_analysis/correl/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/13_spectral_analysis/correl
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./correl > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./correl > omp4_output.txt 2>&1
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
| **Conversion order** | 162 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (14 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `four1`, `realft`, `twofft` |
| **Tags** | `spectral-analysis`, `signal-processing`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
