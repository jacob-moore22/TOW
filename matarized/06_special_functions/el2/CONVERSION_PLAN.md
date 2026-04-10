---
example: el2
chapter: 06_special_functions
chapter_title: "Special Functions"
status: not_started
complexity: medium
conversion_order: 164
priority: 0
tags: [special-function, mathematical-function, cross-chapter]
dependencies: [qsimp, trapzd, ran3]
reverse_dependencies: []
---

# EL2 -- Special Functions

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `el2` |
| **Chapter** | 06 -- Special Functions |
| **Purpose** | General elliptic integral of the second kind. |
| **Status** | `not_started` |
| **Complexity** | `medium` |
| **Fortran LOC** | 50 |
| **Subroutine** | `EL2` (function) |

## 2. Source Files

- **Fortran source:** `fortran/06_special_functions/el2/el2.f` (50 lines)
- **Driver/demo:** `fortran/06_special_functions/el2/el2.dem`
- **Target:** `matarized/06_special_functions/el2/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `qsimp` (04_integration)
  - `trapzd` (04_integration)
  - `ran3` (07_random_numbers)

### Diagram

```mermaid
graph TD
    el2 --> qsimp
    el2 --> trapzd
    el2 --> ran3
```

### Cross-Chapter Dependencies

- `qsimp` from chapter 04
- `trapzd` from chapter 04
- `ran3` from chapter 07

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `AA` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `BB` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `CA` | `REAL` | (scalar) | constant | `constexpr double CA = .0003;` | constant = .0003 |
| `CB` | `REAL` | (scalar) | constant | `constexpr double CB = 1.E-9;` | constant = 1.E-9 |
| `PI` | `REAL` | (scalar) | constant | `constexpr double PI = 3.14159265;` | constant = 3.14159265 |
| `QQC` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X` | `REAL` | (scalar) | parameter (input) | `double` |  |

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
inline double el2(double x, double qqc, double aa, double bb)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `EL2` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(el2_matar_parallel CXX)

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
set(RANDOMNUMBERS_DIR    ${MATARIZED_ROOT}/07_random_numbers)

# --- Build the EL2 example ---
add_executable(el2 main.cpp)
target_link_libraries(el2 matar_lib)
target_include_directories(el2 PRIVATE
    ${INTEGRATION_DIR}/qsimp
    ${INTEGRATION_DIR}/trapzd
    ${RANDOMNUMBERS_DIR}/ran3
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
cd fortran/06_special_functions/el2
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/06_special_functions/el2
mkdir -p build && cd build
cmake .. && make
./el2 > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/06_special_functions/el2/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/06_special_functions/el2
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./el2 > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./el2 > omp4_output.txt 2>&1
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
| **Conversion order** | 164 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | medium (50 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `qsimp`, `trapzd`, `ran3` |
| **Tags** | `special-function`, `mathematical-function`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
