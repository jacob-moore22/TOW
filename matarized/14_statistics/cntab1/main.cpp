#include <cstdio>
#include <cmath>
#include <matar.h>
#include "cntab1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NI = 2, NJ = 2;
        DFMatrixKokkos<double> nn(NI, NJ);

        // Independent 2x2 table: equal proportions
        nn.host(1, 1) = 50.0; nn.host(1, 2) = 50.0;
        nn.host(2, 1) = 50.0; nn.host(2, 2) = 50.0;

        double chisq, df, prob, cramrv, ccc;
        cntab1(nn, NI, NJ, chisq, df, prob, cramrv, ccc);

        std::printf("Contingency table (independent 2x2):\n");
        std::printf("  chisq  = %14.6f  (expected ~ 0)\n", chisq);
        std::printf("  df     = %8.1f\n", df);
        std::printf("  prob   = %14.6f  (expected ~ 1)\n", prob);
        std::printf("  Cramer = %14.6f\n", cramrv);
        std::printf("  CCC    = %14.6f\n", ccc);

        // Strongly associated table
        nn.host(1, 1) = 90.0; nn.host(1, 2) = 10.0;
        nn.host(2, 1) = 10.0; nn.host(2, 2) = 90.0;

        cntab1(nn, NI, NJ, chisq, df, prob, cramrv, ccc);

        std::printf("\nContingency table (strongly associated):\n");
        std::printf("  chisq  = %14.6f\n", chisq);
        std::printf("  df     = %8.1f\n", df);
        std::printf("  prob   = %14.10f\n", prob);
        std::printf("  Cramer = %14.6f\n", cramrv);
        std::printf("  CCC    = %14.6f\n", ccc);

        bool pass = chisq > 50.0 && prob < 0.001;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
