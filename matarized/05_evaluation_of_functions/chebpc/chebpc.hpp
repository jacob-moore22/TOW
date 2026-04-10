#pragma once
#include <matar.h>

using namespace mtr;

// Convert Chebyshev coefficients c(1..n) to polynomial coefficients d(1..n).
// The polynomial is in the shifted variable y = (2x - a - b) / (b - a).
inline void chebpc(const DFMatrixKokkos<double>& c,
                   DFMatrixKokkos<double>& d, int n)
{
    DFMatrixKokkos<double> dd(n);
    for (int j = 1; j <= n; j++)
        dd.host(j) = 0.0;

    for (int j = 1; j <= n; j++)
        d.host(j) = 0.0;

    d.host(1) = c.host(n);

    for (int j = n - 1; j >= 2; j--) {
        for (int k = n - j + 1; k >= 2; k--) {
            double sv = d.host(k);
            d.host(k)  = 2.0 * d.host(k - 1) - dd.host(k);
            dd.host(k) = sv;
        }
        double sv = d.host(1);
        d.host(1)  = -dd.host(1) + c.host(j);
        dd.host(1) = sv;
    }

    for (int j = n; j >= 2; j--)
        d.host(j) = d.host(j - 1) - dd.host(j);

    d.host(1) = -dd.host(1) + 0.5 * c.host(1);
}
