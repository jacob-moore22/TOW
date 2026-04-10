#include <cstdio>
#include <cmath>
#include <matar.h>
#include "tred2.hpp"

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
        a.host(1,1) = 2.0;  a.host(1,2) = 1.0;  a.host(1,3) = 0.0;
        a.host(2,1) = 1.0;  a.host(2,2) = 3.0;  a.host(2,3) = 1.0;
        a.host(3,1) = 0.0;  a.host(3,2) = 1.0;  a.host(3,3) = 2.0;

        tred2(a, N, d, e);

        std::printf("Householder tridiagonal reduction (TRED2)\n");
        std::printf("Input: [[2,1,0],[1,3,1],[0,1,2]]\n\n");

        std::printf("Diagonal d:\n");
        for (int i = 1; i <= N; i++)
            std::printf("  d(%d) = %16.10f\n", i, d.host(i));

        std::printf("\nOff-diagonal e:\n");
        for (int i = 1; i <= N; i++)
            std::printf("  e(%d) = %16.10f\n", i, e.host(i));

        std::printf("\nOrthogonal transformation Q (columns):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.6f", a.host(i, j));
            std::printf(" ]\n");
        }

        // Verify: trace of tridiagonal = trace of original = 2+3+2 = 7
        double trace = 0.0;
        for (int i = 1; i <= N; i++)
            trace += d.host(i);

        std::printf("\nTrace of tridiagonal: %.10f (expected 7.0)\n", trace);
        std::printf("Test %s\n", std::fabs(trace - 7.0) < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
