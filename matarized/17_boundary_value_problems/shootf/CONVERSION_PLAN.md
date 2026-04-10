---
example: shootf
chapter: 17_boundary_value_problems
chapter_title: "Boundary Value Problems"
status: not_started
complexity: high
conversion_order: 186
priority: 0
tags: [boundary-value, differential-equation, cross-chapter]
dependencies: [lubksb, ludcmp, odeint, rk4, rkqc]
reverse_dependencies: []
---

# SHOOTF -- Boundary Value Problems

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `shootf` |
| **Chapter** | 17 -- Boundary Value Problems |
| **Purpose** | Solve a two-point BVP by shooting from both boundaries to a fitting point. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 57 |
| **Subroutine** | `SHOOTF` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/17_boundary_value_problems/shootf/shootf.f` (57 lines)
- **Driver/demo:** `fortran/17_boundary_value_problems/shootf/shootf.dem`
- **Target:** `matarized/17_boundary_value_problems/shootf/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `lubksb` (02_linear_algebra)
  - `ludcmp` (02_linear_algebra)
  - `odeint` (16_ode_integration)
  - `rk4` (16_ode_integration)
  - `rkqc` (16_ode_integration)

### Diagram

```mermaid
graph TD
    shootf --> lubksb
    shootf --> ludcmp
    shootf --> odeint
    shootf --> rk4
    shootf --> rkqc
```

### Cross-Chapter Dependencies

- `lubksb` from chapter 02
- `ludcmp` from chapter 02
- `odeint` from chapter 16
- `rk4` from chapter 16
- `rkqc` from chapter 16

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `DELV1` | `REAL` | N2 | parameter (input) | `DFMatrixKokkos<double>(N2)` |  |
| `DELV2` | `REAL` | N1 | parameter (input) | `DFMatrixKokkos<double>(N1)` |  |
| `DFDV` | `REAL` | NP, NP | local | `DFMatrixKokkos<double>(NP, NP)` |  |
| `DV1` | `REAL` | N2 | parameter (input) | `DFMatrixKokkos<double>(N2)` |  |
| `DV2` | `REAL` | N1 | parameter (input) | `DFMatrixKokkos<double>(N1)` |  |
| `EPS` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `F` | `REAL` | NVAR | parameter (input) | `DFMatrixKokkos<double>(NVAR)` |  |
| `F1` | `REAL` | NP | local | `DFMatrixKokkos<double>(NP)` |  |
| `F2` | `REAL` | NP | local | `DFMatrixKokkos<double>(NP)` |  |
| `H1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `HMI N` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `INDX` | `INTEGER` | NP | local | `DFMatrixKokkos<int>(NP)` |  |
| `N1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `N2` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NP` | `INTEGER` | (scalar) | constant | `constexpr int NP = 20;` | constant = 20 |
| `NVAR` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `V1` | `REAL` | N2 | parameter (input) | `DFMatrixKokkos<double>(N2)` |  |
| `V2` | `REAL` | N1 | parameter (input) | `DFMatrixKokkos<double>(N1)` |  |
| `X1` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `X2` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `XF` | `REAL` | (scalar) | parameter (input) | `double` |  |
| `Y` | `REAL` | NP | local | `DFMatrixKokkos<double>(NP)` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 12  IV=1,N2

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K2: DO 11  I=1,NVAR

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K3: DO 14  IV=1,N1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K4: DO 13  I=1,NVAR

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K5: DO 15  I=1,NVAR

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K6: DO 16  IV=1,N2

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J

### K7: DO 17  IV=1,N1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: J


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
inline void shootf(int nvar, DFMatrixKokkos<double>& v1, DFMatrixKokkos<double>& v2, DFMatrixKokkos<double>& delv1, DFMatrixKokkos<double>& delv2, int n1, int n2, double x1, double x2, double xf, double eps, double h1, double hmi n, DFMatrixKokkos<double>& f, DFMatrixKokkos<double>& dv1, DFMatrixKokkos<double>& dv2)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `SHOOTF` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(shootf_matar_parallel CXX)

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
set(ODEINTEGRATION_DIR   ${MATARIZED_ROOT}/16_ode_integration)

# --- Build the SHOOTF example ---
add_executable(shootf main.cpp)
target_link_libraries(shootf matar_lib)
target_include_directories(shootf PRIVATE
    ${LINEARALGEBRA_DIR}/lubksb
    ${LINEARALGEBRA_DIR}/ludcmp
    ${ODEINTEGRATION_DIR}/odeint
    ${ODEINTEGRATION_DIR}/rk4
    ${ODEINTEGRATION_DIR}/rkqc
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
cd fortran/17_boundary_value_problems/shootf
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/17_boundary_value_problems/shootf
mkdir -p build && cd build
cmake .. && make
./shootf > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/17_boundary_value_problems/shootf/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/17_boundary_value_problems/shootf
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./shootf > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./shootf > omp4_output.txt 2>&1
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
| **Conversion order** | 186 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (57 Fortran LOC, 5 dependencies) |
| **Prerequisite conversions** | `lubksb`, `ludcmp`, `odeint`, `rk4`, `rkqc` |
| **Tags** | `boundary-value`, `differential-equation`, `cross-chapter` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
