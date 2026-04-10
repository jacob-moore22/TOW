#pragma once
#include <vector>
#include <matar.h>

using namespace mtr;

// Convert Chebyshev coefficients c(1..n) to polynomial coefficients d(1..n).
// The polynomial is in the shifted variable y = (2x - a - b) / (b - a).
inline void chebpc(const DFMatrixKokkos<double>& c,
                   DFMatrixKokkos<double>& d, int n)
{
    std::vector<double> dd(n + 1, 0.0);

    for (int j = 1; j <= n; j++)
        d(j) = 0.0;

    d(1) = c(n);

    for (int j = n - 1; j >= 2; j--) {
        for (int k = n - j + 1; k >= 2; k--) {
            double sv = d(k);
            d(k)  = 2.0 * d(k - 1) - dd[k];
            dd[k] = sv;
        }
        double sv = d(1);
        d(1)  = -dd[1] + c(j);
        dd[1] = sv;
    }

    for (int j = n; j >= 2; j--)
        d(j) = d(j - 1) - dd[j];

    d(1) = -dd[1] + 0.5 * c(1);
}
