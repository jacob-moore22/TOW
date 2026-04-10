#pragma once
#include <matar.h>

using namespace mtr;

// Polynomial division u/v. Returns quotient q and remainder r.
// Coefficients stored as: (1)=constant, (2)=x^1, etc.
inline void poldiv(const DFMatrixKokkos<double>& u, int n,
                   const DFMatrixKokkos<double>& v, int nv,
                   DFMatrixKokkos<double>& q,
                   DFMatrixKokkos<double>& r)
{
    for (int j = 1; j <= n; j++) {
        r(j) = u(j);
        q(j) = 0.0;
    }

    for (int k = n - nv; k >= 0; k--) {
        q(k + 1) = r(nv + k) / v(nv);
        for (int j = nv + k - 1; j >= k + 1; j--)
            r(j) -= q(k + 1) * v(j - k);
    }

    r(nv) = 0.0;
}
