#include <cstdio>
#include <cmath>
#include <matar.h>

using namespace mtr;

// Pull in solvde and its helpers
#include "../bksub/bksub.hpp"
#include "../pinvs/pinvs.hpp"
#include "../red/red.hpp"
#include "solvde.hpp"
#include "../difeq/difeq.hpp"

// Test: solve the Legendre equation eigenvalue problem n=2, m=0
// The eigenvalue should be n*(n+1) = 6.
int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NE  = 3;
        constexpr int M   = 41;
        constexpr int NB  = 1;
        constexpr int NCI = NE;
        constexpr int NCJ = NE - NB + 1;
        constexpr int NCK = M + 1;
        constexpr int NSI = NE;
        constexpr int NSJ = 2 * NE + 1;

        int mm_val = 0, n_val = 2;
        double c2 = 0.0;

        double h = 1.0 / (M - 1);
        double x[M + 1];
        for (int k = 1; k <= M - 1; k++)
            x[k] = (k - 1) * h;
        x[M] = 1.0;

        double anorm = 1.0;

        int indexv[4];
        if ((n_val + mm_val) % 2 == 1) {
            indexv[1] = 1; indexv[2] = 2; indexv[3] = 3;
        } else {
            indexv[1] = 2; indexv[2] = 1; indexv[3] = 3;
        }

        DFMatrixKokkos<double> y(NE, M);
        DFMatrixKokkos<double> c(NCI, NCJ, NCK);
        DFMatrixKokkos<double> s(NSI, NSJ);

        for (int k = 1; k <= M; k++) {
            y.host(1, k) = 1.0;
            y.host(2, k) = 0.0;
            y.host(3, k) = static_cast<double>(n_val * (n_val + 1) - mm_val * (mm_val + 1));
        }

        SpheroidalDifeq difeq;
        difeq.x      = x;
        difeq.h      = h;
        difeq.mm     = mm_val;
        difeq.n      = n_val;
        difeq.c2     = c2;
        difeq.anorm  = anorm;
        difeq.m_grid = M;

        double scalv[4];
        scalv[1] = std::fabs(anorm);
        scalv[2] = std::fabs(anorm);
        scalv[3] = std::fmax(1.0, static_cast<double>(n_val * (n_val + 1)));

        std::printf("solvde test: spheroidal harmonics n=%d, m=%d, c2=%.1f\n",
                    n_val, mm_val, c2);

        solvde(100, 5.0e-6, 1.0, scalv, indexv, NE, NB, M,
               y, NE, M, c, NCI, NCJ, NCK, s, NSI, NSJ, difeq);

        double lambda = y.host(3, 1) + mm_val * (mm_val + 1);
        double expected = static_cast<double>(n_val * (n_val + 1));
        double error = std::fabs(lambda - expected);

        std::printf("Lambda = %.10f  (expected %.1f, error = %.2e)\n",
                    lambda, expected, error);
        std::printf("Test %s\n", error < 1e-4 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
