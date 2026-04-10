#include <cstdio>
#include <cmath>
#include <matar.h>
#include "polcof.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Recover coefficients of p(x) = 2 + 3x + x^2
        constexpr int N = 3;
        DFMatrixKokkos<double> xa(N);
        DFMatrixKokkos<double> ya(N);
        DFMatrixKokkos<double> cof(N);

        double xvals[] = {1.0, 2.0, 3.0};
        double exact_cof[] = {2.0, 3.0, 1.0};

        for (int i = 1; i <= N; i++) {
            xa.host(i) = xvals[i - 1];
            double xi  = xa.host(i);
            ya.host(i) = 2.0 + 3.0 * xi + xi * xi;
        }

        polcof(xa, ya, N, cof);

        std::printf("Polynomial coefficients via synthetic division\n");
        std::printf("for p(x) = 2 + 3x + x^2\n\n");

        std::printf("%6s %14s %14s %14s\n", "k", "cof(k)", "exact", "error");
        double max_err = 0.0;
        for (int k = 1; k <= N; k++) {
            double err = std::fabs(cof.host(k) - exact_cof[k - 1]);
            if (err > max_err) max_err = err;
            std::printf("%6d %14.8f %14.8f %14.2e\n",
                        k, cof.host(k), exact_cof[k - 1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
