#include <cstdio>
#include <cmath>
#include <matar.h>
#include "svdcmp.hpp"
#include "svbksb.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int M = 3, N = 3;

        DFMatrixKokkos<double> a(M, N);
        DFMatrixKokkos<double> w(N);
        DFMatrixKokkos<double> v(N, N);
        DFMatrixKokkos<double> b(M);
        DFMatrixKokkos<double> x(N);

        double A_vals[3][3] = {{1,2,3},{4,5,6},{7,8,10}};
        double b_vals[3]    = {14, 32, 53};
        double x_exact[3]   = {1.0, 2.0, 3.0};

        for (int i = 1; i <= M; i++) {
            for (int j = 1; j <= N; j++)
                a.host(i, j) = A_vals[i-1][j-1];
            b.host(i) = b_vals[i-1];
        }

        svdcmp(a, M, N, w, v);
        svbksb(a, w, v, M, N, b, x);

        std::printf("SVD Back-Substitution (SVBKSB)\n");
        std::printf("===============================\n\n");

        std::printf("Solution:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(x.host(i) - x_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  x(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, x.host(i), x_exact[i-1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
