#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "sort.hpp"
#include "probks.hpp"

using namespace mtr;

// Kolmogorov-Smirnov one-sample test. Compares data(1..n) against
// the cumulative distribution function func(x).
// Returns KS statistic d and significance prob.
inline void ksone(DFMatrixKokkos<double>& data, int n,
                  double (*func)(double),
                  double& d, double& prob)
{
    sort(n, data);
    double en = n;
    d = 0.0;
    double fo = 0.0;

    for (int j = 1; j <= n; j++) {
        double fn = j / en;
        double ff = func(data.host(j));
        double dt = std::max(fabs(fo - ff), fabs(fn - ff));
        if (dt > d) d = dt;
        fo = fn;
    }

    prob = probks(sqrt(en) * d);
}
