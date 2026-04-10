// Objective function for robust line fitting (Numerical Recipes ROFUNC).
//
// Given slope b, computes the sum of x(j)*sign(residual(j)).
// Also returns the median intercept aa and the sum of absolute
// deviations abdev. Called by medfit for bisection root-finding.
// Replaces the Fortran COMMON block with explicit parameters.

#pragma once
#include <cmath>
#include <matar.h>
#include "sort.hpp"

using namespace mtr;

inline double rofunc(double b, int ndata,
                     DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                     DFMatrixKokkos<double>& arr,
                     double& aa, double& abdev)
{
    for (int j = 1; j <= ndata; j++)
        arr.host(j) = y.host(j) - b * x.host(j);

    sort(ndata, arr);

    int n1  = ndata + 1;
    int nml = n1 / 2;
    int nmh = n1 - nml;
    aa = 0.5 * (arr.host(nml) + arr.host(nmh));

    double sum = 0.0;
    abdev = 0.0;
    for (int j = 1; j <= ndata; j++) {
        double d = y.host(j) - (b * x.host(j) + aa);
        abdev += std::fabs(d);
        sum += x.host(j) * std::copysign(1.0, d);
    }
    return sum;
}
