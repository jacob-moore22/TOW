#include <cstdio>
#include <cmath>
#include <matar.h>
#include "svdvar.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int MA = 3;
        DFMatrixKokkos<double> v(MA, MA);
        DFMatrixKokkos<double> w(MA);
        DFMatrixKokkos<double> cvm(MA, MA);

        // Identity V matrix with known singular values
        for (int i = 1; i <= MA; i++)
            for (int j = 1; j <= MA; j++)
                v.host(i, j) = (i == j) ? 1.0 : 0.0;

        w.host(1) = 2.0;
        w.host(2) = 3.0;
        w.host(3) = 5.0;

        svdvar(v, MA, w, cvm);

        // With V=I, cvm(i,j) = delta(i,j) / w(i)^2
        std::printf("Covariance Matrix from SVD (SVDVAR)\n");
        std::printf("===================================\n\n");
        std::printf("V = I,  w = {%.1f, %.1f, %.1f}\n\n", w.host(1), w.host(2), w.host(3));

        std::printf("Covariance matrix:\n");
        double max_err = 0.0;
        for (int i = 1; i <= MA; i++) {
            std::printf("  ");
            for (int j = 1; j <= MA; j++) {
                std::printf("%12.6f", cvm.host(i, j));
                double expected = (i == j) ? 1.0 / (w.host(i) * w.host(i)) : 0.0;
                double err = std::fabs(cvm.host(i, j) - expected);
                if (err > max_err) max_err = err;
            }
            std::printf("\n");
        }

        std::printf("\nExpected diagonal: {%.6f, %.6f, %.6f}\n",
                    1.0/4.0, 1.0/9.0, 1.0/25.0);
        std::printf("Max error: %.2e\n", max_err);

        // Verify symmetry
        double sym_err = 0.0;
        for (int i = 1; i <= MA; i++)
            for (int j = 1; j <= MA; j++)
                sym_err = std::fmax(sym_err, std::fabs(cvm.host(i, j) - cvm.host(j, i)));

        std::printf("Symmetry error: %.2e\n", sym_err);
        std::printf("Test %s\n", (max_err < 1e-12 && sym_err < 1e-12) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
