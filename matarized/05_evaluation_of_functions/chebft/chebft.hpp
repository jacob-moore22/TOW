#pragma once
#include <cmath>
#include <vector>
#include <matar.h>

using namespace mtr;

// Chebyshev fit: compute coefficients c(1..n) approximating func on [a,b].
template<typename Func>
inline void chebft(double a, double b, DFMatrixKokkos<double>& c, int n, Func func)
{
    constexpr double PI = 3.141592653589793;

    std::vector<double> f(n + 1);

    double bma = 0.5 * (b - a);
    double bpa = 0.5 * (b + a);

    for (int k = 1; k <= n; k++) {
        double y = std::cos(PI * (k - 0.5) / n);
        f[k] = func(y * bma + bpa);
    }

    double fac = 2.0 / n;
    for (int j = 1; j <= n; j++) {
        double sum = 0.0;
        for (int k = 1; k <= n; k++) {
            sum += f[k] * std::cos((PI * (j - 1)) * ((k - 0.5) / n));
        }
        c(j) = fac * sum;
    }
}
