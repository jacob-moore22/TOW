// Sparse matrix iterative solver using conjugate-gradient on normal
// equations (Numerical Recipes SPARSE).
//
// Solves A*x = b where A is only available through matrix-vector products.
//
// b(1..n)  -- right-hand side
// n        -- system size
// asub     -- function computing asub(x, out) = A*x
// atsub    -- function computing atsub(x, out) = A^T*x
// x(1..n)  -- initial guess on input, solution on output
// rsq      -- residual squared on output
//
// asub/atsub are std::function taking (DFMatrixKokkos<double>&, DFMatrixKokkos<double>&).

#pragma once
#include <cmath>
#include <cstdio>
#include <functional>
#include <matar.h>

using namespace mtr;

using SparseMV = std::function<void(DFMatrixKokkos<double>&, DFMatrixKokkos<double>&)>;

inline void sparse(DFMatrixKokkos<double>& b, int n,
                   SparseMV asub, SparseMV atsub,
                   DFMatrixKokkos<double>& x, double& rsq)
{
    constexpr double EPS = 1.0e-6;
    double eps2 = n * EPS * EPS;

    DFMatrixKokkos<double> g(n), h(n), xi(n), xj(n);

    int irst = 0;
restart:
    irst++;
    asub(x, xi);

    double rp = 0.0, bsq = 0.0;
    for (int j = 1; j <= n; j++) {
        bsq += b.host(j) * b.host(j);
        xi.host(j) -= b.host(j);
        rp += xi.host(j) * xi.host(j);
    }

    atsub(xi, g);
    for (int j = 1; j <= n; j++) {
        g.host(j) = -g.host(j);
        h.host(j) =  g.host(j);
    }

    for (int iter = 1; iter <= 10 * n; iter++) {
        asub(h, xi);

        double anum = 0.0, aden = 0.0;
        for (int j = 1; j <= n; j++) {
            anum += g.host(j) * h.host(j);
            aden += xi.host(j) * xi.host(j);
        }
        if (aden == 0.0) {
            std::fprintf(stderr, "sparse: very singular matrix\n");
            return;
        }
        anum /= aden;

        for (int j = 1; j <= n; j++) {
            xi.host(j) = x.host(j);
            x.host(j) += anum * h.host(j);
        }

        asub(x, xj);

        rsq = 0.0;
        for (int j = 1; j <= n; j++) {
            xj.host(j) -= b.host(j);
            rsq += xj.host(j) * xj.host(j);
        }

        if (rsq == rp || rsq <= bsq * eps2) return;

        if (rsq > rp) {
            for (int j = 1; j <= n; j++)
                x.host(j) = xi.host(j);
            if (irst >= 3) return;
            goto restart;
        }

        rp = rsq;
        atsub(xj, xi);

        double gg = 0.0, dgg = 0.0;
        for (int j = 1; j <= n; j++) {
            gg  += g.host(j) * g.host(j);
            dgg += (xi.host(j) + g.host(j)) * xi.host(j);
        }
        if (gg == 0.0) return;

        double gam = dgg / gg;
        for (int j = 1; j <= n; j++) {
            g.host(j) = -xi.host(j);
            h.host(j) =  g.host(j) + gam * h.host(j);
        }
    }

    std::fprintf(stderr, "sparse: too many iterations\n");
}
