#include <cstdio>
#include <cmath>
#include <matar.h>
#include "mmid.hpp"

using namespace mtr;

void harmonic_derivs(double /*x*/, DFMatrixKokkos<double>& y,
                     DFMatrixKokkos<double>& dydx)
{
    dydx.host(1) =  y.host(2);
    dydx.host(2) = -y.host(1);
}

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NVAR = 2;
        DFMatrixKokkos<double> y(NVAR), dydx(NVAR), yout(NVAR);

        double xs   = 0.0;
        double htot = 1.0;

        y.host(1) = 0.0;
        y.host(2) = 1.0;
        harmonic_derivs(xs, y, dydx);

        std::printf("MMID convergence test: y'' = -y, step from 0 to %.1f\n\n",
                    htot);
        std::printf("%8s %16s %16s %12s %12s\n",
                    "nstep", "y1", "sin(1)", "err_y1", "err_y2");

        int steps[] = {2, 4, 6, 8, 12, 16, 24, 48};
        double max_err = 1.0;
        for (int s = 0; s < 8; s++) {
            mmid(y, dydx, NVAR, xs, htot, steps[s], yout, harmonic_derivs);
            double err1 = std::fabs(yout.host(1) - std::sin(htot));
            double err2 = std::fabs(yout.host(2) - std::cos(htot));
            max_err = std::max(err1, err2);
            std::printf("%8d %16.10f %16.10f %12.2e %12.2e\n",
                        steps[s], yout.host(1), std::sin(htot), err1, err2);
        }

        std::printf("\nBest result (nstep=48):\n");
        std::printf("  y1 = %16.10f  (exact sin(1) = %16.10f)\n",
                    yout.host(1), std::sin(1.0));
        std::printf("  y2 = %16.10f  (exact cos(1) = %16.10f)\n",
                    yout.host(2), std::cos(1.0));
        std::printf("  Max error: %.2e\n", max_err);
        std::printf("  Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
