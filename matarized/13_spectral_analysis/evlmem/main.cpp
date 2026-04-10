// Driver for MATAR evlmem -- evaluate MEM power spectrum.
// Uses known AR(1) coefficients and verifies the analytic spectrum.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "evlmem.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // AR(1) process: x(n) = a*x(n-1) + e,  power = var_e / |1 - a*exp(-i*w)|^2
        const int M = 1;
        const double a = 0.9;
        const double pm = 0.19;  // ~ 1 - a^2

        DFMatrixKokkos<double> cof(M);
        cof.host(1) = a;
        cof.update_device();

        std::printf("MEM spectrum of AR(1) with a=%.2f, pm=%.4f\n", a, pm);
        std::printf("  fdt       EVLMEM      Analytic\n");

        double max_err = 0.0;
        for (int i = 0; i <= 20; i++) {
            double fdt = 0.5 * i / 20.0;
            double result = evlmem(fdt, cof, M, pm);

            double theta = 6.28318530717959 * fdt;
            double sr = 1.0 - a * std::cos(theta);
            double si = -a * std::sin(theta);
            double analytic = pm / (sr * sr + si * si);

            double err = std::fabs(result - analytic);
            if (err > max_err) max_err = err;

            std::printf("  %.3f   %12.6f %12.6f\n", fdt, result, analytic);
        }
        std::printf("\n  Max |evlmem - analytic| = %.2e\n\n", max_err);
    }
    MATAR_FINALIZE();
    return 0;
}
