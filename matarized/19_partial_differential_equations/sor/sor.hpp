#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Successive over-relaxation for solving elliptic PDEs
// (Numerical Recipes SOR).
//
// Solves: a*u(j+1,l) + b*u(j-1,l) + c*u(j,l+1) + d*u(j,l-1) + e*u(j,l) = f(j,l)
// on interior points j,l = 2..jmax-1.
// Red-black ordering enables parallelism via FOR_ALL.
// rjac: spectral radius of the Jacobi iteration (typically cos(pi/jmax)).
inline void sor(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b,
                DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& d,
                DFMatrixKokkos<double>& e, DFMatrixKokkos<double>& f,
                DFMatrixKokkos<double>& u, int jmax, double rjac)
{
    const int    MAXITS = 10000;
    const double EPS    = 1.0e-5;

    double loc_anormf = 0.0, anormf = 0.0;
    DO_REDUCE_SUM(j, 2, jmax - 1,
                  l, 2, jmax - 1,
                  loc_anormf, {
        loc_anormf += fabs(f(j, l));
    }, anormf);

    double omega = 1.0;

    for (int n = 1; n <= MAXITS; n++) {
        int parity = n % 2;

        double loc_anorm = 0.0, anorm = 0.0;
        DO_REDUCE_SUM(j, 2, jmax - 1,
                      l, 2, jmax - 1,
                      loc_anorm, {
            if ((j + l) % 2 == parity) {
                double resid = a(j, l) * u(j + 1, l)
                             + b(j, l) * u(j - 1, l)
                             + c(j, l) * u(j, l + 1)
                             + d(j, l) * u(j, l - 1)
                             + e(j, l) * u(j, l)
                             - f(j, l);
                loc_anorm += fabs(resid);
                u(j, l)    = u(j, l) - omega * resid / e(j, l);
            }
        }, anorm);

        if (n == 1)
            omega = 1.0 / (1.0 - 0.5 * rjac * rjac);
        else
            omega = 1.0 / (1.0 - 0.25 * rjac * rjac * omega);

        if (n > 1 && anorm < EPS * anormf) {
            std::printf("  SOR converged in %d iterations\n", n);
            return;
        }
    }
    std::fprintf(stderr, "sor: MAXITS exceeded\n");
}
