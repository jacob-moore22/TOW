#include <cstdio>
#include <cmath>
#include <matar.h>
#include "rkdumb.hpp"

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
        constexpr int NVAR  = 2;
        constexpr int NSTEP = 200;
        DFMatrixKokkos<double> vstart(NVAR);
        DFMatrixKokkos<double> xx(NSTEP + 1);
        DFMatrixKokkos<double> y_out(NVAR, NSTEP + 1);

        vstart.host(1) = 0.0;
        vstart.host(2) = 1.0;

        rkdumb(vstart, NVAR, 0.0, 2.0 * M_PI, NSTEP, harmonic_derivs,
               xx, y_out);

        std::printf("RKDUMB integration of y'' = -y from 0 to 2*pi (%d steps)\n\n",
                    NSTEP);
        std::printf("%10s %16s %16s %16s %16s\n",
                    "x", "y1", "sin(x)", "y2", "cos(x)");

        for (int k = 1; k <= NSTEP + 1; k += NSTEP / 4) {
            double xk = xx.host(k);
            std::printf("%10.4f %16.10f %16.10f %16.10f %16.10f\n",
                        xk, y_out.host(1, k), std::sin(xk),
                        y_out.host(2, k), std::cos(xk));
        }

        int last = NSTEP + 1;
        double err1 = std::fabs(y_out.host(1, last));
        double err2 = std::fabs(y_out.host(2, last) - 1.0);
        double max_err = std::max(err1, err2);

        std::printf("\ny1(2*pi) = %16.10f  (exact 0)\n", y_out.host(1, last));
        std::printf("y2(2*pi) = %16.10f  (exact 1)\n", y_out.host(2, last));
        std::printf("Max error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
