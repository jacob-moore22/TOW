---
example: twofft
chapter: 12_fourier_transform
chapter_title: "Fourier Transform"
status: converted
complexity: low
conversion_order: 54
priority: 2
tags: [fourier-transform, fft, signal-processing]
dependencies: [four1]
reverse_dependencies: [convlv, correl]
---

# TWOFFT -- Fourier Transform

## 1. Overview

| Field | Value |
|-------|-------|
| **Example** | `twofft` |
| **Chapter** | 12 -- Fourier Transform |
| **Purpose** | Simultaneous FFT of two real functions packed into one complex FFT. |
| **Status** | `converted` |
| **Complexity** | `low` |
| **Fortran LOC** | 23 |
| **Subroutine** | `TWOFFT` (subroutine) |

## 2. Source Files

- **Fortran source:** `fortran/12_fourier_transform/twofft/twofft.f` (23 lines)
- **Driver/demo:** `fortran/12_fourier_transform/twofft/twofft.dem`
- **Target:** `matarized/12_fourier_transform/twofft/`


## 3. Dependency Graph

### Forward Dependencies (this example depends on)

  - `four1` (12_fourier_transform)

### Diagram

```mermaid
graph TD
    twofft --> four1
```

### Cross-Chapter Dependencies

(none)

## 4. Reverse Dependencies (examples that depend on this)

  - `convlv` (13_spectral_analysis)
  - `correl` (13_spectral_analysis)

> **Conversion note:** This routine is depended on by 2 other examples and should be converted early.

## 5. Fortran Variable Catalog

| Name | Fortran Type | Shape | Role | MATAR Type | Notes |
|------|-------------|-------|------|-----------|-------|
| `C1` | `COMPLEX` | (scalar) | local | `double` |  |
| `C2` | `COMPLEX` | (scalar) | local | `double` |  |
| `DATA1` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `DATA2` | `REAL` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `FFT1` | `COMPLEX` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `FFT2` | `COMPLEX` | N | parameter (input) | `DFMatrixKokkos<double>(N)` |  |
| `H1` | `COMPLEX` | (scalar) | local | `double` |  |
| `H2` | `COMPLEX` | (scalar) | local | `double` |  |
| `N` | `INTEGER` | (scalar) | parameter (input) | `int` |  |

### MATAR Type Mapping Rationale

- **Layout:** `FMatrix` (column-major) preserves Fortran memory layout for correctness.
- **Index base:** `Matrix` (1-based) matches Fortran indexing with `DO_ALL` inclusive ranges.
- **Residence:** `Dual` (`DFMatrixKokkos`) enables both host I/O and device computation.
- **Ownership:** Owning types at call site; consider `ViewFMatrix` for sub-array slices.

## 6. Compute Kernel Analysis

### K1: DO 11  J=1,N  (interleave two real arrays into one complex array)

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Each iteration writes `DATA1(2*J-1)`, `DATA1(2*J)` at positions uniquely determined by J. Embarrassingly parallel.

### K2: DO 12  J=2,N/2+1  (separate the two transforms)

- **Thread safety:** `safe`
- **Recommended macro:** `DO_ALL`
- **Notes:** Each iteration computes FFT1 and FFT2 at unique complex indices `(J)` and `(N+2-J)`. The two writes per iteration are at distinct addresses. Embarrassingly parallel.


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
inline void twofft(DFMatrixKokkos<double>& data1, DFMatrixKokkos<double>& data2, DFMatrixKokkos<double>& fft1, DFMatrixKokkos<double>& fft2, int n)
```

### Output Format

- **.hpp header** (included by other examples via `#include`)

### Steps

1. **Translate data structures** -- replace Fortran arrays with `DFMatrixKokkos` (see variable catalog below)
2. **Translate routine** -- convert `TWOFFT` to a C++ function as a `.hpp header`
3. **Replace loops** -- convert DO loops to `DO_ALL` / `DO_REDUCE_*` macros (see kernel analysis below)
4. **Add synchronization** -- insert `MATAR_FENCE()` between dependent kernels; add `update_host()`/`update_device()` for Dual types
5. **Create driver** -- translate the `.dem` test program to `main.cpp` with `MATAR_INITIALIZE` / `MATAR_FINALIZE` boilerplate
6. **Generate CMakeLists.txt** -- use the template below (based on convlv reference)
7. **Validate** -- follow the validation plan below

## 8. CMake Configuration

Based on the [convlv CMakeLists.txt](../../13_spectral_analysis/convlv/CMakeLists.txt) reference template.

```cmake
cmake_minimum_required(VERSION 3.18)
project(twofft_matar_parallel CXX)

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

# --- Build the TWOFFT example ---
add_executable(twofft main.cpp)
target_link_libraries(twofft matar_lib)
target_include_directories(twofft PRIVATE
    ${FOURIERTRANSFORM_DIR}/four1
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
cd fortran/12_fourier_transform/twofft
make run > reference_output.txt 2>&1
```


### Serial Validation

```bash
cd matarized/12_fourier_transform/twofft
mkdir -p build && cd build
cmake .. && make
./twofft > serial_output.txt 2>&1
diff <(head -50 serial_output.txt) <(head -50 ../../../../fortran/12_fourier_transform/twofft/reference_output.txt)
```


### Parallel Validation (OpenMP)

```bash
cd matarized/12_fourier_transform/twofft
mkdir -p build-omp && cd build-omp
cmake .. -DENABLE_OPENMP=ON && make
OMP_NUM_THREADS=1 ./twofft > omp1_output.txt 2>&1
OMP_NUM_THREADS=4 ./twofft > omp4_output.txt 2>&1
# Verify: omp1 output must exactly match serial output
diff serial_output.txt omp1_output.txt
# Verify: omp4 output must match within floating-point tolerance
```


### Pass Criteria

- Max absolute difference vs. Fortran reference: **< 1e-10** (double precision)

- OpenMP results must be deterministic across repeated runs

- No runtime errors, memory leaks, or Kokkos warnings


## Conversion Audit (Already Converted)

This example has been converted to MATAR at `matarized/12_fourier_transform/twofft/twofft.hpp`.

### Key Design Decisions

- **Interleave/separate loops** are both embarrassingly parallel with `DO_ALL`.
- **Calls four1:** Depends on `four1.hpp` for the underlying complex FFT.
- **No sequential dependencies:** Unlike four1 and realft, twofft has no twiddle recurrences.

### Remaining Work

- [ ] Create standalone `main.cpp` driver with test cases from `twofft.dem`
- [ ] Add CMakeLists.txt for independent builds

## 11. Agent Metadata

| Field | Value |
|-------|-------|
| **Conversion order** | 54 of 202 |
| **Priority score** | 2 (reverse dependency count) |
| **Estimated effort** | low (23 Fortran LOC, 1 dependencies) |
| **Prerequisite conversions** | `four1` |
| **Tags** | `fourier-transform`, `fft`, `signal-processing` |
| **MATAR reference sections** | Sec 5 (parallel loops), Sec 15 (Fortran interop) |
