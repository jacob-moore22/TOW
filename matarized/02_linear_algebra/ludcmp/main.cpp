#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ludcmp.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<int>    indx(N);

        // Same 3x3 system as gaussj test
        double A_vals[3][3] = {{1,2,3},{4,5,6},{7,8,10}};

        for (int i = 1; i <= N; i++)
            for (int j = 1; j <= N; j++)
                a.host(i, j) = A_vals[i-1][j-1];

        double d;
        ludcmp(a, N, indx, d);

        std::printf("LU Decomposition (LUDCMP)\n");
        std::printf("=========================\n\n");

        std::printf("LU matrix (combined L and U):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  ");
            for (int j = 1; j <= N; j++)
                std::printf("%12.6f", a.host(i, j));
            std::printf("\n");
        }

        std::printf("\nPivot indices: ");
        for (int i = 1; i <= N; i++)
            std::printf("%d ", indx.host(i));
        std::printf("\n");

        std::printf("Determinant sign: %.0f\n", d);

        // Compute determinant = d * product of diagonal
        double det = d;
        for (int i = 1; i <= N; i++)
            det *= a.host(i, i);
        std::printf("Determinant: %.6f\n", det);

        // For the matrix [[1,2,3],[4,5,6],[7,8,10]], det = -3
        double det_exact = -3.0;
        double err = std::fabs(det - det_exact);
        std::printf("Expected determinant: %.1f\n", det_exact);
        std::printf("Determinant error: %.2e\n", err);
        std::printf("Test %s\n", err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
