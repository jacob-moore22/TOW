#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "tqli.hpp"
#include "../tred2/tred2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<double> d(N);
        DFMatrixKokkos<double> e(N);

        // Symmetric test matrix: [[2,1,0],[1,3,1],[0,1,2]]
        // Known eigenvalues: 1, 2, 4
        a.host(1,1) = 2.0;  a.host(1,2) = 1.0;  a.host(1,3) = 0.0;
        a.host(2,1) = 1.0;  a.host(2,2) = 3.0;  a.host(2,3) = 1.0;
        a.host(3,1) = 0.0;  a.host(3,2) = 1.0;  a.host(3,3) = 2.0;

        std::printf("TRED2 + TQLI eigenvalue decomposition\n");
        std::printf("Input: [[2,1,0],[1,3,1],[0,1,2]]\n");
        std::printf("Expected eigenvalues: 1, 2, 4\n\n");

        tred2(a, N, d, e);

        std::printf("After TRED2:\n");
        std::printf("  Diagonal d:     ");
        for (int i = 1; i <= N; i++) std::printf(" %10.6f", d.host(i));
        std::printf("\n  Off-diagonal e: ");
        for (int i = 1; i <= N; i++) std::printf(" %10.6f", e.host(i));
        std::printf("\n\n");

        tqli(d, e, N, a);

        std::printf("After TQLI:\n");
        std::printf("Computed eigenvalues:\n");
        double evals[N];
        for (int i = 1; i <= N; i++) {
            evals[i - 1] = d.host(i);
            std::printf("  d(%d) = %16.10f\n", i, d.host(i));
        }

        std::printf("\nEigenvectors (columns of Z):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.6f", a.host(i, j));
            std::printf(" ]\n");
        }

        std::sort(evals, evals + N);
        double expected[] = {1.0, 2.0, 4.0};
        double max_err = 0.0;
        for (int i = 0; i < N; i++) {
            double err = std::fabs(evals[i] - expected[i]);
            if (err > max_err) max_err = err;
        }

        std::printf("\nMax eigenvalue error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
