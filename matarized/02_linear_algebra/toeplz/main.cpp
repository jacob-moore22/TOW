#include <cstdio>
#include <cmath>
#include <matar.h>
#include "toeplz.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Symmetric Toeplitz system of size 4:
        // T = [ 4  3  2  1 ]   x = [1]   y = T*x = [10]
        //     [ 3  4  3  2 ]       [1]              [12]
        //     [ 2  3  4  3 ]       [1]              [12]
        //     [ 1  2  3  4 ]       [1]              [10]
        constexpr int N = 4;

        // r has 2*N-1 = 7 elements: r(1..7), with r(N)=r(4) on diagonal
        // Row of Toeplitz: [1, 2, 3, 4, 3, 2, 1]
        DFMatrixKokkos<double> r(2 * N - 1);
        DFMatrixKokkos<double> x(N);
        DFMatrixKokkos<double> y(N);

        double r_vals[] = {1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0};
        double y_vals[] = {10.0, 12.0, 12.0, 10.0};
        double x_exact[] = {1.0, 1.0, 1.0, 1.0};

        for (int i = 1; i <= 2 * N - 1; i++)
            r.host(i) = r_vals[i - 1];
        for (int i = 1; i <= N; i++)
            y.host(i) = y_vals[i - 1];

        toeplz(r, x, y, N);

        std::printf("Toeplitz Solver (TOEPLZ)\n");
        std::printf("========================\n\n");

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
