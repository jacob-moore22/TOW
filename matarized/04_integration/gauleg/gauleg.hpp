#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Compute abscissas x[0..n-1] and weights w[0..n-1] for
// n-point Gauss-Legendre quadrature over [x1, x2].
// Output arrays must be pre-allocated DCArrayKokkos of length n.
inline void gauleg(double x1, double x2,
                   DCArrayKokkos<double>& x, DCArrayKokkos<double>& w, int n)
{
    constexpr double eps = 3.0e-14;
    constexpr double pi  = 3.14159265358979323846;

    int m = (n + 1) / 2;
    double xm = 0.5 * (x2 + x1);
    double xl = 0.5 * (x2 - x1);

    for (int i = 0; i < m; i++) {
        double z = std::cos(pi * (i + 0.75) / (n + 0.5));
        double pp, z1;

        do {
            double p1 = 1.0;
            double p2 = 0.0;
            for (int j = 1; j <= n; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z  = z1 - p1 / pp;
        } while (std::fabs(z - z1) > eps);

        x.host(i)         = xm - xl * z;
        x.host(n - 1 - i) = xm + xl * z;
        w.host(i)         = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        w.host(n - 1 - i) = w.host(i);
    }

    x.update_device();
    w.update_device();
}
