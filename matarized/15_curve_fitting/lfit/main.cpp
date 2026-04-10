#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gaussj.hpp"
#include "covsrt.hpp"
#include "fpoly.hpp"
#include "lfit.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 50;
        constexpr int MA    = 3;
        constexpr int MFIT  = 3;

        DFMatrixKokkos<double> x(NDATA), y(NDATA), sig(NDATA);
        DFMatrixKokkos<double> a(MA);
        DFMatrixKokkos<int>    lista(MA);
        DFMatrixKokkos<double> covar(MA, MA);
        double chisq;

        // Generate y = 1 + 2*x + 3*x^2 with noise, x in [-1, 1]
        unsigned seed = 54321;
        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = -1.0 + 2.0 * (i - 1) / (NDATA - 1);
            seed = seed * 1103515245 + 12345;
            double noise = 0.05 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(i) = 1.0 + 2.0 * x.host(i) + 3.0 * x.host(i) * x.host(i) + noise;
            sig.host(i) = 0.1;
        }

        for (int j = 1; j <= MA; j++) {
            lista.host(j) = j;
            a.host(j) = 0.0;
        }

        auto funcs = [](double xv, DFMatrixKokkos<double>& p, int np) {
            fpoly(xv, p, np);
        };

        lfit(x, y, sig, NDATA, a, MA, lista, MFIT, covar, chisq, funcs);

        std::printf("General Linear Least Squares (LFIT)\n");
        std::printf("===================================\n\n");
        std::printf("Fitting y = a1 + a2*x + a3*x^2  (true: 1 + 2x + 3x^2)\n\n");

        double true_a[] = {1.0, 2.0, 3.0};
        double max_err = 0.0;
        for (int j = 1; j <= MA; j++) {
            double err = std::fabs(a.host(j) - true_a[j - 1]);
            if (err > max_err) max_err = err;
            std::printf("  a(%d) = %10.6f  (true = %.1f, err = %.2e)\n",
                        j, a.host(j), true_a[j - 1], err);
        }

        std::printf("\nChi-squared = %.4f\n", chisq);

        std::printf("\nCovariance matrix:\n");
        for (int i = 1; i <= MA; i++) {
            std::printf("  ");
            for (int j = 1; j <= MA; j++)
                std::printf("%12.6f", covar.host(i, j));
            std::printf("\n");
        }

        std::printf("\nMax coefficient error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 0.1 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
