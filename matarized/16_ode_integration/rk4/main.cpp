#include <cstdio>
#include <cmath>
#include <matar.h>
#include "rk4.hpp"

using namespace mtr;

void harmonic_derivs(double /*x*/, DFMatrixKokkos<double>& y,
                     DFMatrixKokkos<double>& dydx)
{
    dydx.host(1) =  y.host(2);   // y1' =  y2
    dydx.host(2) = -y.host(1);   // y2' = -y1
}

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NVAR = 2;
        DFMatrixKokkos<double> y(NVAR), dydx(NVAR);

        y.host(1) = 0.0;
        y.host(2) = 1.0;

        constexpr int NSTEP = 100;
        double h = 2.0 * M_PI / NSTEP;
        double x = 0.0;

        std::printf("RK4 integration of y'' = -y from 0 to 2*pi (%d steps, h=%.4f)\n\n",
                    NSTEP, h);
        std::printf("%10s %16s %16s %16s %16s\n",
                    "x", "y1", "sin(x)", "y2", "cos(x)");

        for (int k = 1; k <= NSTEP; k++) {
            harmonic_derivs(x, y, dydx);
            rk4(y, dydx, NVAR, x, h, y, harmonic_derivs);
            x = k * h;
            if (k % 25 == 0 || k == NSTEP)
                std::printf("%10.4f %16.10f %16.10f %16.10f %16.10f\n",
                            x, y.host(1), std::sin(x), y.host(2), std::cos(x));
        }

        double err1 = std::fabs(y.host(1) - std::sin(2.0 * M_PI));
        double err2 = std::fabs(y.host(2) - std::cos(2.0 * M_PI));
        double max_err = std::max(err1, err2);

        std::printf("\ny1(2*pi) = %16.10f  (exact 0)\n", y.host(1));
        std::printf("y2(2*pi) = %16.10f  (exact 1)\n", y.host(2));
        std::printf("Max error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
