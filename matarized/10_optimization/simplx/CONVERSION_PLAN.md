---
example: simplx
chapter: 10_optimization
chapter_title: "Optimization"
status: not_started
complexity: high
conversion_order: 171
priority: 0
tags: [optimization, minimization]
dependencies: [simp1, simp2, simp3]
reverse_dependencies: []
---

# SIMPLX -- Optimization

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `simplx` |
| **Chapter** | 10 -- Optimization |
| **Purpose** | Linear programming using the simplex method. |
| **Status** | `not_started` |
| **Complexity** | `high` |
| **Fortran LOC** | 96 |
| **Subroutine** | `SIMPLX` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/10_optimization/simplx/simplx.f` (96 lines)
- **Driver/demo:** `fortran/10_optimization/simplx/simplx.dem`
- **Target:** `matarized/10_optimization/simplx/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `simp1` (10_optimization)
  - `simp2` (10_optimization)
  - `simp3` (10_optimization)

### Diagram

```mermaid
graph TD
    simplx --> simp1
    simplx --> simp2
    simplx --> simp3
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  (none)

> **Conversion note:** No other examples depend on this routine.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `A` | `REAL` | MP, NP | parameter (input) | `DFMatrixKokkos<double>(MP, NP)` |  |
| `EPS` | `REAL` | (scalar) | constant | `constexpr double EPS = 1.E-6;` | constant = 1.E-6 |
| `ICASE` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `IPOSV` | `INTEGER` | M | parameter (input) | `DFMatrixKokkos<int>(M)` |  |
| `IZROV` | `INTEGER` | N | parameter (input) | `DFMatrixKokkos<int>(N)` |  |
| `L1` | `INTEGER` | MMAX | local | `DFMatrixKokkos<int>(MMAX)` |  |
| `L2` | `INTEGER` | MMAX | local | `DFMatrixKokkos<int>(MMAX)` |  |
| `L3` | `INTEGER` | MMAX | local | `DFMatrixKokkos<int>(MMAX)` |  |
| `M` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `M1` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `M2` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `M3` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `MMAX` | `INTEGER` | (scalar) | constant | `constexpr int MMAX = 100;` | constant = 100 |
| `MP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |
| `NP` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  K=1,N

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K2: DO 12  I=1,M

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K3: DO 13  I=1,M2

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K4: DO 15  K=1,N+1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K5: DO 14  I=M1+1,M

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K6: DO 16  IP=M12,M

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K7: DO 18  I=M1+1,M12

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K8: DO 17  K=1,N+1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K9: DO 19  K=1,NL1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K10: DO 21  IS=K,NL1

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1

### K11: DO 22  I=1,M+2

- **Thread safety:** `reduction`
- **Recommended macro:** `DO_REDUCE_SUM`
- **Notes:** Accumulates: Q1


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
inline void simplx(DFMatrixKokkos<double>& a, int m, int n, int mp, int np, int m1, int m2, int m3, int icase, DFMatrixKokkos<double>& izrov, DFMatrixKokkos<double>& iposv)
```

### Output Format

- **.cpp with main()** (standalone executable)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `SIMPLX` to a C++ function as a `.cpp with main()`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(simplx_matar_parallel CXX)

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
set(OPTIMIZATION_DIR     ${MATARIZED_ROOT}/10_optimization)

# --- Build the SIMPLX example ---
add_executable(simplx main.cpp)
target_link_libraries(simplx matar_lib)
target_include_directories(simplx PRIVATE
    ${OPTIMIZATION_DIR}/simp1
    ${OPTIMIZATION_DIR}/simp2
    ${OPTIMIZATION_DIR}/simp3
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
cd fortran/10_optimization/simplx
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/10_optimization/simplx
mkdir -p build && cd build
cmake .. && make
./simplx > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/10_optimization/simplx/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/10_optimization/simplx
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./simplx > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./simplx > omp4_output.txt 2>&1
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
| **Conversion order** | 171 of 202 |
| **Priority score** | 0 (reverse dependency count) |
| **Estimated effort** | high (96 Fortran LOC, 3 dependencies) |
| **Prerequisite conversions** | `simp1`, `simp2`, `simp3` |
| **Tags** | `optimization`, `minimization` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 6 (reductions), Sec 15 (Fortran interop) |
