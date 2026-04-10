#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gaussj.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;
        constexpr int M = 1;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<double> b(N, M);

        // 3x3 system:  A * x = b
        //  [ 1  2  3 ] [x1]   [14]
        //  [ 4  5  6 ] [x2] = [32]
        //  [ 7  8 10 ] [x3]   [53]
        // Known solution: x = [1, 2, 3]
        double A_vals[3][3] = {{1,2,3},{4,5,6},{7,8,10}};
        double b_vals[3]    = {14, 32, 53};
        double x_exact[3]   = {1.0, 2.0, 3.0};

        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++)
                a.host(i, j) = A_vals[i-1][j-1];
            b.host(i, 1) = b_vals[i-1];
        }

        // Save original A for residual check
        double A_orig[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                A_orig[i][j] = A_vals[i][j];

        gaussj(a, N, b, M);

        std::printf("Gauss-Jordan Elimination (GAUSSJ)\n");
        std::printf("==================================\n\n");

        std::printf("Solution:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(b.host(i, 1) - x_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  x(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, b.host(i, 1), x_exact[i-1], err);
        }

        std::printf("\nInverse matrix (A^-1):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  ");
            for (int j = 1; j <= N; j++)
                std::printf("%12.6f", a.host(i, j));
            std::printf("\n");
        }

        // Verify A^-1 * A_orig = I
        std::printf("\nVerification A^-1 * A_orig = I:\n");
        double id_err = 0.0;
        for (int i = 0; i < N; i++) {
            std::printf("  ");
            for (int j = 0; j < N; j++) {
                double sum = 0.0;
                for (int k = 0; k < N; k++)
                    sum += a.host(i+1, k+1) * A_orig[k][j];
                double expected = (i == j) ? 1.0 : 0.0;
                id_err = std::fmax(id_err, std::fabs(sum - expected));
                std::printf("%12.6f", sum);
            }
            std::printf("\n");
        }

        std::printf("\nMax solution error:  %.2e\n", max_err);
        std::printf("Max identity error:  %.2e\n", id_err);
        std::printf("Test %s\n", (max_err < 1e-10 && id_err < 1e-10) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
