#include <cstdio>
#include <cmath>
#include <matar.h>
#include "laguer.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== LAGUER: Laguerre's method for polynomial roots ===\n\n");
        bool all_pass = true;

        // p(x) = x^2 - 2 → roots at +/- sqrt(2)
        {
            const int m = 2;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1) = -2.0; a.host(2) = 0.0;   // coef 0: -2 + 0i
            a.host(3) =  0.0; a.host(4) = 0.0;   // coef 1:  0 + 0i
            a.host(5) =  1.0; a.host(6) = 0.0;   // coef 2:  1 + 0i

            double xr = 1.0, xi = 0.0;
            laguer(a, m, xr, xi, 1.0e-12, true);

            double err = std::fabs(xr - std::sqrt(2.0)) + std::fabs(xi);
            std::printf("x^2 - 2:  root = (%.15f, %.15f)\n", xr, xi);
            std::printf("  expected sqrt(2) = %.15f, err = %.2e\n",
                        std::sqrt(2.0), err);
            if (err > 1e-10) all_pass = false;
        }

        // p(x) = x^3 - x - 1 → one real root near 1.3247
        {
            const int m = 3;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1) = -1.0; a.host(2) = 0.0;
            a.host(3) = -1.0; a.host(4) = 0.0;
            a.host(5) =  0.0; a.host(6) = 0.0;
            a.host(7) =  1.0; a.host(8) = 0.0;

            double xr = 1.5, xi = 0.0;
            laguer(a, m, xr, xi, 1.0e-12, true);

            double val = xr * xr * xr - xr - 1.0;
            std::printf("\nx^3-x-1:  root = (%.15f, %.15f)\n", xr, xi);
            std::printf("  f(root) = %.2e\n", std::fabs(val));
            if (std::fabs(val) > 1e-10) all_pass = false;
        }

        // p(x) = x^2 + 1 → roots at +/- i
        {
            const int m = 2;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1) = 1.0; a.host(2) = 0.0;
            a.host(3) = 0.0; a.host(4) = 0.0;
            a.host(5) = 1.0; a.host(6) = 0.0;

            double xr = 0.5, xi = 0.5;
            laguer(a, m, xr, xi, 1.0e-12, true);

            double err = std::fabs(xr) + std::fabs(std::fabs(xi) - 1.0);
            std::printf("\nx^2 + 1:  root = (%.15f, %.15f)\n", xr, xi);
            std::printf("  expected (0, +/-1), err = %.2e\n", err);
            if (err > 1e-10) all_pass = false;
        }

        std::printf("\nTest %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
