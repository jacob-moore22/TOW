#include <cstdio>
#include <cmath>
#include <matar.h>
#include "balanc.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);

        // Asymmetric matrix with large scaling differences
        a.host(1,1) = 1.0;    a.host(1,2) = 100.0;  a.host(1,3) = 0.01;
        a.host(2,1) = 0.001;  a.host(2,2) = 2.0;    a.host(2,3) = 50.0;
        a.host(3,1) = 10.0;   a.host(3,2) = 0.005;  a.host(3,3) = 3.0;

        std::printf("Matrix balancing (BALANC)\n\n");
        std::printf("Before balancing:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %12.6f", a.host(i, j));
            std::printf(" ]\n");
        }

        double norm_before = 0.0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                norm_before += a.host(i, j) * a.host(i, j);
        norm_before = std::sqrt(norm_before);

        balanc(a, N);

        std::printf("\nAfter balancing:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %12.6f", a.host(i, j));
            std::printf(" ]\n");
        }

        double norm_after = 0.0;
        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                norm_after += a.host(i, j) * a.host(i, j);
        norm_after = std::sqrt(norm_after);

        // Row/column norms should be more balanced
        std::printf("\nRow norms after balancing:\n");
        for (int i = 1; i <= N; i++) {
            double rn = 0.0;
            for (int j = 1; j <= N; j++)
                rn += std::fabs(a.host(i, j));
            std::printf("  Row %d: %12.6f\n", i, rn);
        }

        std::printf("\nColumn norms after balancing:\n");
        for (int i = 1; i <= N; i++) {
            double cn = 0.0;
            for (int j = 1; j <= N; j++)
                cn += std::fabs(a.host(j, i));
            std::printf("  Col %d: %12.6f\n", i, cn);
        }

        std::printf("\nFrobenius norm before: %12.6f\n", norm_before);
        std::printf("Frobenius norm after:  %12.6f\n", norm_after);
        std::printf("Test %s\n", norm_after <= norm_before ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
