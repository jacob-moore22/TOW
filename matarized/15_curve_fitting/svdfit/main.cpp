#include <cstdio>
#include <cmath>
#include <matar.h>
#include "svdcmp.hpp"
#include "svbksb.hpp"
#include "fpoly.hpp"
#include "svdfit.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 50;
        constexpr int MA    = 3;

        DFMatrixKokkos<double> x(NDATA), y(NDATA), sig(NDATA);
        DFMatrixKokkos<double> a(MA);
        DFMatrixKokkos<double> u(NDATA, MA), v(MA, MA), w(MA);
        double chisq;

        // Generate y = 1 + 2*x + 3*x^2 with noise, x in [-1, 1]
        unsigned seed = 67890;
        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = -1.0 + 2.0 * (i - 1) / (NDATA - 1);
            seed = seed * 1103515245 + 12345;
            double noise = 0.05 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(i) = 1.0 + 2.0 * x.host(i) + 3.0 * x.host(i) * x.host(i) + noise;
            sig.host(i) = 0.1;
        }

        auto funcs = [](double xv, DFMatrixKokkos<double>& p, int np) {
            fpoly(xv, p, np);
        };

        svdfit(x, y, sig, NDATA, a, MA, u, v, w, chisq, funcs);

        std::printf("SVD Least Squares Fit (SVDFIT)\n");
        std::printf("==============================\n\n");
        std::printf("Fitting y = a1 + a2*x + a3*x^2  (true: 1 + 2x + 3x^2)\n\n");

        double true_a[] = {1.0, 2.0, 3.0};
        double max_err = 0.0;
        for (int j = 1; j <= MA; j++) {
            double err = std::fabs(a.host(j) - true_a[j - 1]);
            if (err > max_err) max_err = err;
            std::printf("  a(%d) = %10.6f  (true = %.1f, err = %.2e)\n",
                        j, a.host(j), true_a[j - 1], err);
        }

        std::printf("\nSingular values: ");
        for (int j = 1; j <= MA; j++)
            std::printf("%.6f  ", w.host(j));
        std::printf("\n");

        std::printf("Chi-squared = %.4f\n", chisq);
        std::printf("\nMax coefficient error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 0.1 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
