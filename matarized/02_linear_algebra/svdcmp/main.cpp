#include <cstdio>
#include <cmath>
#include <matar.h>
#include "svdcmp.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int M = 3, N = 3;

        DFMatrixKokkos<double> a(M, N);
        DFMatrixKokkos<double> w(N);
        DFMatrixKokkos<double> v(N, N);

        double A_orig[3][3] = {{1,2,3},{4,5,6},{7,8,10}};

        for (int i = 1; i <= M; i++)
            for (int j = 1; j <= N; j++)
                a.host(i, j) = A_orig[i-1][j-1];

        svdcmp(a, M, N, w, v);

        std::printf("Singular Value Decomposition (SVDCMP)\n");
        std::printf("======================================\n\n");

        std::printf("Singular values:\n");
        for (int i = 1; i <= N; i++)
            std::printf("  w(%d) = %16.10f\n", i, w.host(i));

        std::printf("\nU matrix:\n");
        for (int i = 1; i <= M; i++) {
            std::printf("  ");
            for (int j = 1; j <= N; j++)
                std::printf("%12.6f", a.host(i, j));
            std::printf("\n");
        }

        std::printf("\nV matrix:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  ");
            for (int j = 1; j <= N; j++)
                std::printf("%12.6f", v.host(i, j));
            std::printf("\n");
        }

        // Verify: reconstruct A = U * diag(W) * V^T
        std::printf("\nReconstruction A = U*W*V^T:\n");
        double max_err = 0.0;
        for (int i = 1; i <= M; i++) {
            std::printf("  ");
            for (int j = 1; j <= N; j++) {
                double sum = 0.0;
                for (int k = 1; k <= N; k++)
                    sum += a.host(i, k) * w.host(k) * v.host(j, k);
                double err = std::fabs(sum - A_orig[i-1][j-1]);
                if (err > max_err) max_err = err;
                std::printf("%12.6f", sum);
            }
            std::printf("\n");
        }

        std::printf("\nOriginal matrix:\n");
        for (int i = 0; i < M; i++) {
            std::printf("  ");
            for (int j = 0; j < N; j++)
                std::printf("%12.6f", A_orig[i][j]);
            std::printf("\n");
        }

        std::printf("\nMax reconstruction error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
