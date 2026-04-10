#include <cstdio>
#include <cmath>
#include <matar.h>
#include "difeq.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NE = 3;
        constexpr int M  = 11;
        constexpr int NSI = NE;
        constexpr int NSJ = 2 * NE + 1;

        DFMatrixKokkos<double> s(NSI, NSJ);
        DFMatrixKokkos<double> y(NE, M);

        double x[M + 1];
        double h = 1.0 / (M - 1);
        for (int k = 1; k <= M; k++)
            x[k] = (k - 1) * h;

        for (int j = 1; j <= NE; j++)
            for (int k = 1; k <= M; k++)
                y.host(j, k) = 1.0;

        int indexv[4] = {0, 1, 2, 3};

        SpheroidalDifeq difeq;
        difeq.x      = x;
        difeq.h      = h;
        difeq.mm     = 0;
        difeq.n      = 2;
        difeq.c2     = 0.0;
        difeq.anorm  = 1.0;
        difeq.m_grid = M;

        std::printf("difeq test: interior point k=5\n");
        difeq(5, 1, M, NSJ, 1, NE, indexv, NE, s, NSI, NSJ, y, NE, M);

        std::printf("s matrix after difeq(k=5):\n");
        for (int i = 1; i <= NSI; i++) {
            for (int j = 1; j <= NSJ; j++)
                std::printf(" %8.4f", s.host(i, j));
            std::printf("\n");
        }

        std::printf("\ndifeq test: boundary k=1\n");
        difeq(1, 1, M, NSJ, 1, NE, indexv, NE, s, NSI, NSJ, y, NE, M);

        std::printf("s matrix after difeq(k=1):\n");
        for (int i = 1; i <= NSI; i++) {
            for (int j = 1; j <= NSJ; j++)
                std::printf(" %8.4f", s.host(i, j));
            std::printf("\n");
        }

        std::printf("Test PASSED (difeq completed without error)\n");
    }
    MATAR_FINALIZE();
    return 0;
}
