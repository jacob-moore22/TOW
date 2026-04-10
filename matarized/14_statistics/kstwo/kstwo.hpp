#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "sort.hpp"
#include "probks.hpp"

using namespace mtr;

// Kolmogorov-Smirnov two-sample test. Compares data1(1..n1) and data2(1..n2).
// Returns KS statistic d and significance prob.
inline void kstwo(DFMatrixKokkos<double>& data1, int n1,
                  DFMatrixKokkos<double>& data2, int n2,
                  double& d, double& prob)
{
    sort(n1, data1);
    sort(n2, data2);

    double en1 = n1;
    double en2 = n2;
    int j1 = 1, j2 = 1;
    double fo1 = 0.0, fo2 = 0.0;
    d = 0.0;

    while (j1 <= n1 && j2 <= n2) {
        if (data1.host(j1) < data2.host(j2)) {
            double fn1 = j1 / en1;
            double dt = std::max(fabs(fn1 - fo2), fabs(fo1 - fo2));
            if (dt > d) d = dt;
            fo1 = fn1;
            j1++;
        } else {
            double fn2 = j2 / en2;
            double dt = std::max(fabs(fn2 - fo1), fabs(fo2 - fo1));
            if (dt > d) d = dt;
            fo2 = fn2;
            j2++;
        }
    }

    prob = probks(sqrt(en1 * en2 / (en1 + en2)) * d);
}
