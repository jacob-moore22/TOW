---
example: flmoon
chapter: 01_dates_and_calendars
chapter_title: "Dates and Calendars"
status: not_started
complexity: medium
conversion_order: 92
priority: 1
tags: [utility, calendar]
dependencies: [caldat, julday]
reverse_dependencies: [badluk]
---

# FLMOON -- Dates and Calendars

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `flmoon` |
| **Chapter** | 01 -- Dates and Calendars |
| **Purpose** | Calculate date and time of phases of the moon. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 26 |
| **Subroutine** | `FLMOON` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/01_dates_and_calendars/flmoon/flmoon.f` (26 lines)
- **Driver/demo:** `fortran/01_dates_and_calendars/flmoon/flmoon.dem`
- **Target:** `matarized/01_dates_and_calendars/flmoon/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `caldat` (01_dates_and_calendars)
  - `julday` (01_dates_and_calendars)

### Diagram

```mermaid
graph TD
    flmoon --> caldat
    flmoon --> julday
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `badluk` (01_dates_and_calendars)

> **Conversion note:** This routine is depended on by 1 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `FRAC` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `JD` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NPH` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `RAD` | `REAL` | (scalar) | constant | `constexpr double RAD = 0.017453293;` | constant = 0.017453293 |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

_No DO loops detected in source._

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
inline void flmoon(int n, int nph, int jd, double frac)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `FLMOON` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(flmoon_matar_parallel CXX)

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
set(DATESANDCALENDARS_DIR ${MATARIZED_ROOT}/01_dates_and_calendars)

# --- Build the FLMOON example ---
add_executable(flmoon main.cpp)
target_link_libraries(flmoon matar_lib)
target_include_directories(flmoon PRIVATE
    ${DATESANDCALENDARS_DIR}/caldat
    ${DATESANDCALENDARS_DIR}/julday
)
```

## 9. Performance Improvements

- **FMatrix to CArray migration:** The initial translation uses `DFMatrixKokkos` (column-major, 1-based) for Fortran compatibility.  For GPU targets, converting to `DCArrayKokkos` (row-major, 0-based) with reordered loops will improve coalesced memory access.
- **Reduction fusion:** If multiple reductions share the same loop bounds, consider fusing them into a single pass to reduce kernel launch overhead.
- **Fence elimination:** After conversion, audit `MATAR_FENCE()` placement.  Remove fences between independent kernels that do not share data.

## 10. Validation Plan

### Reference Output

Build and run the Fortran version to capture reference output:

```bash
cd fortran/01_dates_and_calendars/flmoon
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/01_dates_and_calendars/flmoon
mkdir -p build && cd build
cmake .. && make
./flmoon > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/01_dates_and_calendars/flmoon/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/01_dates_and_calendars/flmoon
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./flmoon > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./flmoon > omp4_output.txt 2>&1
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
| **Conversion order** | 92 of 202 |
| **Priority score** | 1 (reverse dependency count) |
| **Estimated effort** | medium (26 Fortran LOC, 2 dependencies) |
| **Prerequisite conversions** | `caldat`, `julday` |
| **Tags** | `utility`, `calendar` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
