#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gammln.hpp"
#include "gser.hpp"
#include "gcf.hpp"
#include "gammq.hpp"
#include "fit.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 20;
        DFMatrixKokkos<double> x(NDATA);
        DFMatrixKokkos<double> y(NDATA);
        DFMatrixKokkos<double> sig(NDATA);

        // Generate y = 2 + 3*x with deterministic noise
        unsigned seed = 12345;
        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = static_cast<double>(i);
            seed = seed * 1103515245 + 12345;
            double noise = 0.5 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(i) = 2.0 + 3.0 * x.host(i) + noise;
            sig.host(i) = 1.0;
        }

        double a, b, siga, sigb, chi2, q;

        std::printf("Linear Least Squares Fit (FIT)\n");
        std::printf("==============================\n\n");
        std::printf("Fitting y = a + b*x to data from y = 2 + 3*x + noise\n\n");

        // Unweighted fit
        fit(x, y, NDATA, sig, 0, a, b, siga, sigb, chi2, q);

        std::printf("Unweighted fit (mwt=0):\n");
        std::printf("  a    = %10.6f  +/- %8.6f  (true = 2.0)\n", a, siga);
        std::printf("  b    = %10.6f  +/- %8.6f  (true = 3.0)\n", b, sigb);
        std::printf("  chi2 = %10.4f\n", chi2);
        std::printf("  q    = %10.6f\n", q);

        double err_a = std::fabs(a - 2.0);
        double err_b = std::fabs(b - 3.0);

        // Weighted fit
        fit(x, y, NDATA, sig, 1, a, b, siga, sigb, chi2, q);

        std::printf("\nWeighted fit (mwt=1, sig=1):\n");
        std::printf("  a    = %10.6f  +/- %8.6f  (true = 2.0)\n", a, siga);
        std::printf("  b    = %10.6f  +/- %8.6f  (true = 3.0)\n", b, sigb);
        std::printf("  chi2 = %10.4f\n", chi2);
        std::printf("  q    = %10.6f\n", q);

        double err_a2 = std::fabs(a - 2.0);
        double err_b2 = std::fabs(b - 3.0);

        bool pass = (err_a < 1.0) && (err_b < 0.1) && (err_a2 < 1.0) && (err_b2 < 0.1);
        std::printf("\nRecovery: |a-2| = %.4f, |b-3| = %.4f\n", err_a, err_b);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
