#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ludcmp.hpp"
#include "lubksb.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<int>    indx(N);
        DFMatrixKokkos<double> b(N);

        // Same 3x3 system: solution x = [1, 2, 3]
        double A_vals[3][3] = {{1,2,3},{4,5,6},{7,8,10}};
        double b_vals[3]    = {14, 32, 53};
        double x_exact[3]   = {1.0, 2.0, 3.0};

        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++)
                a.host(i, j) = A_vals[i-1][j-1];
            b.host(i) = b_vals[i-1];
        }

        double d;
        ludcmp(a, N, indx, d);
        lubksb(a, N, indx, b);

        std::printf("LU Back-Substitution (LUBKSB)\n");
        std::printf("==============================\n\n");

        std::printf("Solution:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(b.host(i) - x_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  x(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, b.host(i), x_exact[i-1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
