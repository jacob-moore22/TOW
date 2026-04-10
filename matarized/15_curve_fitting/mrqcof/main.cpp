#include <cstdio>
#include <cmath>
#include <matar.h>
#include "fgauss.hpp"
#include "mrqcof.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 20;
        constexpr int MA    = 3;
        constexpr int MFIT  = 3;

        DFMatrixKokkos<double> x(NDATA), y(NDATA), sig(NDATA);
        DFMatrixKokkos<double> a(MA);
        DFMatrixKokkos<int>    lista(MA);
        DFMatrixKokkos<double> alpha(MFIT, MFIT);
        DFMatrixKokkos<double> beta(MFIT);
        double chisq;

        // True Gaussian: A=5, center=10, width=3
        a.host(1) = 5.0;
        a.host(2) = 10.0;
        a.host(3) = 3.0;

        for (int j = 1; j <= MA; j++) lista.host(j) = j;

        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = static_cast<double>(i);
            double arg = (x.host(i) - 10.0) / 3.0;
            y.host(i) = 5.0 * std::exp(-arg * arg);
            sig.host(i) = 0.1;
        }

        auto funcs = [](double xv, DFMatrixKokkos<double>& av, double& yv,
                        DFMatrixKokkos<double>& dyda, int na) {
            fgauss(xv, av, yv, dyda, na);
        };

        mrqcof(x, y, sig, NDATA, a, MA, lista, MFIT, alpha, beta, chisq, funcs);

        std::printf("Normal Equations for LM (MRQCOF)\n");
        std::printf("================================\n\n");

        std::printf("Alpha matrix (curvature):\n");
        for (int i = 1; i <= MFIT; i++) {
            std::printf("  ");
            for (int j = 1; j <= MFIT; j++)
                std::printf("%14.4f", alpha.host(i, j));
            std::printf("\n");
        }

        std::printf("\nBeta vector (gradient):\n");
        for (int j = 1; j <= MFIT; j++)
            std::printf("  beta(%d) = %14.6f\n", j, beta.host(j));

        std::printf("\nChi-squared = %.6f\n", chisq);

        // At the true parameters, beta should be ~0 and alpha should be positive definite
        double beta_norm = 0.0;
        for (int j = 1; j <= MFIT; j++)
            beta_norm += beta.host(j) * beta.host(j);
        beta_norm = std::sqrt(beta_norm);

        bool pos_diag = true;
        for (int j = 1; j <= MFIT; j++)
            if (alpha.host(j, j) <= 0.0) pos_diag = false;

        // Verify alpha symmetry
        double sym_err = 0.0;
        for (int i = 1; i <= MFIT; i++)
            for (int j = 1; j <= MFIT; j++)
                sym_err = std::fmax(sym_err, std::fabs(alpha.host(i, j) - alpha.host(j, i)));

        std::printf("\n|beta| = %.6f (expected ~0 at true params)\n", beta_norm);
        std::printf("Alpha positive diagonal: %s\n", pos_diag ? "yes" : "no");
        std::printf("Alpha symmetry error: %.2e\n", sym_err);
        std::printf("Test %s\n", (beta_norm < 1.0 && pos_diag && sym_err < 1e-10) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
