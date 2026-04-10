#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <matar.h>
#include "sor.hpp"

using namespace mtr;

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
        const double PI = 3.141592653589793;

        if (JMAX < 3) {
            std::fprintf(stderr, "JMAX must be >= 3\n");
            MATAR_FINALIZE();
            return 1;
        }

        DFMatrixKokkos<double> a(JMAX, JMAX), b(JMAX, JMAX), c(JMAX, JMAX),
                               d(JMAX, JMAX), e(JMAX, JMAX), f(JMAX, JMAX),
                               u(JMAX, JMAX);

        // Laplace equation with unit point source at center
        for (int i = 1; i <= JMAX; i++) {
            for (int j = 1; j <= JMAX; j++) {
                a.host(i, j) =  1.0;
                b.host(i, j) =  1.0;
                c.host(i, j) =  1.0;
                d.host(i, j) =  1.0;
                e.host(i, j) = -4.0;
                f.host(i, j) =  0.0;
                u.host(i, j) =  0.0;
            }
        }

        int midl = JMAX / 2 + 1;
        f.host(midl, midl) = 2.0;

        a.update_device();
        b.update_device();
        c.update_device();
        d.update_device();
        e.update_device();
        f.update_device();
        u.update_device();

        double rjac = std::cos(PI / JMAX);

        std::printf("SOR test: Laplace equation on %d x %d grid\n", JMAX, JMAX);
        auto t0 = std::chrono::high_resolution_clock::now();
        sor(a, b, c, d, e, f, u, JMAX, rjac);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        std::printf("  Solve time: %.6f s\n", elapsed);

        MATAR_FENCE();
        u.update_host();

        std::printf("  Center value u(%d,%d) = %.10f\n", midl, midl, u.host(midl, midl));

        // Verify residual
        f.update_host();
        double max_err = 0.0;
        for (int i = 2; i <= JMAX - 1; i++) {
            for (int j = 2; j <= JMAX - 1; j++) {
                double val = u.host(i + 1, j) + u.host(i - 1, j)
                           + u.host(i, j + 1) + u.host(i, j - 1)
                           - 4.0 * u.host(i, j) - f.host(i, j);
                if (std::fabs(val) > std::fabs(max_err))
                    max_err = val;
            }
        }
        std::printf("  Max residual: %.2e\n", max_err);
        std::printf("Test %s\n", std::fabs(max_err) < 1e-4 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
