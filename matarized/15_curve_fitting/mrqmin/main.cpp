#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gaussj.hpp"
#include "covsrt.hpp"
#include "fgauss.hpp"
#include "mrqcof.hpp"
#include "mrqmin.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 100;
        constexpr int MA    = 3;
        constexpr int MFIT  = 3;

        DFMatrixKokkos<double> x(NDATA), y(NDATA), sig(NDATA);
        DFMatrixKokkos<double> a(MA);
        DFMatrixKokkos<int>    lista(MA);
        DFMatrixKokkos<double> covar(MA, MA);
        DFMatrixKokkos<double> alpha(MA, MA);
        double chisq, alamda;

        // True parameters: A=5, center=10, width=3
        double true_amp    = 5.0;
        double true_center = 10.0;
        double true_width  = 3.0;

        unsigned seed = 99999;
        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = 0.2 * i;
            double arg = (x.host(i) - true_center) / true_width;
            seed = seed * 1103515245 + 12345;
            double noise = 0.1 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(i) = true_amp * std::exp(-arg * arg) + noise;
            sig.host(i) = 0.1;
        }

        // Starting guess (perturbed from truth)
        a.host(1) = 4.0;
        a.host(2) = 9.0;
        a.host(3) = 2.5;
        for (int j = 1; j <= MA; j++) lista.host(j) = j;

        auto funcs = [](double xv, DFMatrixKokkos<double>& av, double& yv,
                        DFMatrixKokkos<double>& dyda, int na) {
            fgauss(xv, av, yv, dyda, na);
        };

        std::printf("Levenberg-Marquardt Fit (MRQMIN)\n");
        std::printf("================================\n\n");
        std::printf("True: A=%.1f, center=%.1f, width=%.1f\n", true_amp, true_center, true_width);
        std::printf("Guess: A=%.1f, center=%.1f, width=%.1f\n\n", a.host(1), a.host(2), a.host(3));

        // Initialize
        alamda = -1.0;
        mrqmin(x, y, sig, NDATA, a, MA, lista, MFIT, covar, alpha, chisq, funcs, alamda);
        std::printf("Iter  0: chi2 = %12.4f  alamda = %.6f\n", chisq, alamda);

        // Iterate
        double old_chisq = chisq;
        for (int iter = 1; iter <= 50; iter++) {
            mrqmin(x, y, sig, NDATA, a, MA, lista, MFIT, covar, alpha, chisq, funcs, alamda);
            if (iter <= 10 || iter % 10 == 0)
                std::printf("Iter %2d: chi2 = %12.4f  alamda = %.6f\n", iter, chisq, alamda);
            if (std::fabs(old_chisq - chisq) < 1e-6 * chisq && alamda < 1.0)
                break;
            old_chisq = chisq;
        }

        // Final call for covariance
        alamda = 0.0;
        mrqmin(x, y, sig, NDATA, a, MA, lista, MFIT, covar, alpha, chisq, funcs, alamda);

        std::printf("\nFitted parameters:\n");
        std::printf("  A      = %10.6f  (true = %.1f)\n", a.host(1), true_amp);
        std::printf("  center = %10.6f  (true = %.1f)\n", a.host(2), true_center);
        std::printf("  width  = %10.6f  (true = %.1f)\n", a.host(3), true_width);
        std::printf("  chi2   = %10.4f\n", chisq);

        std::printf("\nCovariance matrix:\n");
        for (int i = 1; i <= MA; i++) {
            std::printf("  ");
            for (int j = 1; j <= MA; j++)
                std::printf("%12.6f", covar.host(i, j));
            std::printf("\n");
        }

        double err_amp    = std::fabs(a.host(1) - true_amp);
        double err_center = std::fabs(a.host(2) - true_center);
        double err_width  = std::fabs(a.host(3) - true_width);

        std::printf("\nRecovery errors:\n");
        std::printf("  |A-5|     = %.4f\n", err_amp);
        std::printf("  |c-10|    = %.4f\n", err_center);
        std::printf("  |w-3|     = %.4f\n", err_width);

        bool pass = (err_amp < 0.5) && (err_center < 0.5) && (err_width < 0.5);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
