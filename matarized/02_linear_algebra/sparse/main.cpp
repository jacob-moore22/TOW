#include <cstdio>
#include <cmath>
#include <matar.h>
#include "sparse.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Solve A*x = b where A is a 4x4 diagonally dominant SPD matrix.
        // A = [ 4 -1  0  0 ]   x_exact = [1]   b = A*x = [ 3]
        //     [-1  4 -1  0 ]              [2]              [ 5]
        //     [ 0 -1  4 -1 ]              [3]              [ 7]
        //     [ 0  0 -1  4 ]              [4]              [13]
        constexpr int N = 4;

        // Store A as dense on host for the matvec callbacks
        double A_vals[4][4] = {
            { 4, -1,  0,  0},
            {-1,  4, -1,  0},
            { 0, -1,  4, -1},
            { 0,  0, -1,  4}
        };

        auto asub = [&](DFMatrixKokkos<double>& xin, DFMatrixKokkos<double>& out) {
            for (int i = 1; i <= N; i++) {
                double s = 0.0;
                for (int j = 1; j <= N; j++)
                    s += A_vals[i-1][j-1] * xin.host(j);
                out.host(i) = s;
            }
        };

        // A is symmetric, so A^T = A
        auto atsub = asub;

        DFMatrixKokkos<double> b(N), x(N);

        double b_vals[]  = {3.0, 5.0, 7.0, 13.0};
        double x_exact[] = {1.0, 2.0, 3.0, 4.0};

        for (int i = 1; i <= N; i++) {
            b.host(i) = b_vals[i-1];
            x.host(i) = 0.0;
        }

        double rsq;
        sparse(b, N, asub, atsub, x, rsq);

        std::printf("Sparse Iterative Solver (SPARSE)\n");
        std::printf("=================================\n\n");

        std::printf("Solution:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(x.host(i) - x_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  x(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, x.host(i), x_exact[i-1], err);
        }

        std::printf("\nResidual squared: %.2e\n", rsq);
        std::printf("Max error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
