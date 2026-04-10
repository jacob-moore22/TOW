#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <matar.h>

using namespace mtr;

#include "../../02_linear_algebra/tridag/tridag.hpp"
#include "adi.hpp"

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int JMAX = 11;
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] != '-') {
                JMAX = std::atoi(argv[i]);
                break;
            }
        }

        if (JMAX < 3) {
            std::fprintf(stderr, "JMAX must be >= 3\n");
            MATAR_FINALIZE();
            return 1;
        }

        // Laplace equation: a=1, b+e = -4 (b=-2, e=-2), c=1, d=1, f=1
        // Splitting: a*u(j-1) + b*u(j) + c*u(j+1)  [x-direction]
        //            d*u(l-1) + e*u(l) + f*u(l+1)  [y-direction]
        DFMatrixKokkos<double> a(JMAX, JMAX), b(JMAX, JMAX), c(JMAX, JMAX),
                               d(JMAX, JMAX), e(JMAX, JMAX), f(JMAX, JMAX),
                               g(JMAX, JMAX), u(JMAX, JMAX);

        for (int i = 1; i <= JMAX; i++) {
            for (int j = 1; j <= JMAX; j++) {
                a.host(i, j) =  1.0;
                b.host(i, j) = -2.0;
                c.host(i, j) =  1.0;
                d.host(i, j) =  1.0;
                e.host(i, j) = -2.0;
                f.host(i, j) =  1.0;
                g.host(i, j) =  0.0;
                u.host(i, j) =  0.0;
            }
        }

        int midl = JMAX / 2 + 1;
        g.host(midl, midl) = -2.0;

        // ADI parameters: k levels, eigenvalue bounds
        int k_adi = 3;
        double alpha = 0.01;
        double beta  = 4.0;
        double eps   = 1.0e-5;

        std::printf("ADI test: Laplace equation on %d x %d grid\n", JMAX, JMAX);
        adi(a, b, c, d, e, f, g, u, JMAX, k_adi, alpha, beta, eps);

        std::printf("  Center value u(%d,%d) = %.10f\n", midl, midl, u.host(midl, midl));

        // Verify residual
        double max_err = 0.0;
        for (int i = 2; i <= JMAX - 1; i++) {
            for (int j = 2; j <= JMAX - 1; j++) {
                double resid = a.host(i, j) * u.host(i - 1, j)
                             + (b.host(i, j) + e.host(i, j)) * u.host(i, j)
                             + c.host(i, j) * u.host(i + 1, j)
                             + d.host(i, j) * u.host(i, j - 1)
                             + f.host(i, j) * u.host(i, j + 1)
                             + g.host(i, j);
                if (std::fabs(resid) > std::fabs(max_err))
                    max_err = resid;
            }
        }
        std::printf("  Max residual: %.2e\n", max_err);
        std::printf("Test %s\n", std::fabs(max_err) < 1e-3 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
