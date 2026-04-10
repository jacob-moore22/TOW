#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bsstep.hpp"

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
        constexpr double EPS  = 1.0e-10;
        constexpr double TINY = 1.0e-30;

        DFMatrixKokkos<double> y(NVAR), dydx(NVAR), yscal(NVAR);

        y.host(1) = 0.0;
        y.host(2) = 1.0;
        double x = 0.0;
        double h = 1.0;

        std::printf("BSSTEP adaptive Bulirsch-Stoer steps for y'' = -y\n\n");
        std::printf("%6s %12s %12s %16s %16s %12s\n",
                    "step", "hdid", "hnext", "y1", "sin(x)", "err");

        for (int k = 1; k <= 12; k++) {
            harmonic_derivs(x, y, dydx);
            for (int i = 1; i <= NVAR; i++)
                yscal.host(i) = std::fabs(y.host(i))
                              + std::fabs(h * dydx.host(i)) + TINY;

            if (x + h > 2.0 * M_PI) h = 2.0 * M_PI - x;

            double hdid, hnext;
            bsstep(y, dydx, NVAR, x, h, EPS, yscal, hdid, hnext,
                   harmonic_derivs);

            double err = std::fabs(y.host(1) - std::sin(x));
            std::printf("%6d %12.6f %12.6f %16.10f %16.10f %12.2e\n",
                        k, hdid, hnext, y.host(1), std::sin(x), err);
            h = hnext;
            if (x >= 2.0 * M_PI - 1.0e-10) break;
        }

        double err1 = std::fabs(y.host(1) - std::sin(x));
        double err2 = std::fabs(y.host(2) - std::cos(x));
        double max_err = std::max(err1, err2);

        std::printf("\nAt x = %.10f:\n", x);
        std::printf("  y1 = %16.10f  (exact sin = %16.10f)\n",
                    y.host(1), std::sin(x));
        std::printf("  y2 = %16.10f  (exact cos = %16.10f)\n",
                    y.host(2), std::cos(x));
        std::printf("  Max error: %.2e\n", max_err);
        std::printf("  Test %s\n", max_err < 1e-8 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
