#include <cstdio>
#include <cmath>
#include <matar.h>
#include "../rkqc/rkqc.hpp"
#include "odeint.hpp"

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
        constexpr int    NVAR  = 2;
        constexpr int    KMAX  = 200;
        constexpr double EPS   = 1.0e-8;
        constexpr double H1    = 0.1;
        constexpr double HMIN  = 1.0e-12;
        constexpr double DXSAV = 0.1;

        DFMatrixKokkos<double> ystart(NVAR);
        ystart.host(1) = 0.0;
        ystart.host(2) = 1.0;

        OdeIntPath path;
        path.kmax  = KMAX;
        path.dxsav = DXSAV;
        path.xp    = DFMatrixKokkos<double>(KMAX);
        path.yp    = DFMatrixKokkos<double>(NVAR, KMAX);

        int nok = 0, nbad = 0;

        auto rkqc_stepper = [](DFMatrixKokkos<double>& y_s,
                               DFMatrixKokkos<double>& dydx_s,
                               int n, double& x_s, double htry, double eps_s,
                               DFMatrixKokkos<double>& yscal_s,
                               double& hdid, double& hnext,
                               auto derivs_fn) {
            rkqc(y_s, dydx_s, n, x_s, htry, eps_s, yscal_s, hdid, hnext,
                 derivs_fn);
        };

        odeint(ystart, NVAR, 0.0, 2.0 * M_PI, EPS, H1, HMIN,
               nok, nbad, harmonic_derivs, rkqc_stepper, path);

        std::printf("ODEINT integration of y'' = -y from 0 to 2*pi\n");
        std::printf("  Stepper: rkqc (adaptive RK4), eps = %.0e\n\n", EPS);
        std::printf("  Steps: %d ok, %d bad, %d stored points\n\n",
                    nok, nbad, path.kount);

        std::printf("%10s %16s %16s %16s %16s\n",
                    "x", "y1", "sin(x)", "y2", "cos(x)");
        int stride = std::max(1, path.kount / 8);
        for (int k = 1; k <= path.kount; k += stride) {
            double xk = path.xp.host(k);
            std::printf("%10.4f %16.10f %16.10f %16.10f %16.10f\n",
                        xk, path.yp.host(1, k), std::sin(xk),
                        path.yp.host(2, k), std::cos(xk));
        }
        if (path.kount % stride != 1) {
            double xk = path.xp.host(path.kount);
            std::printf("%10.4f %16.10f %16.10f %16.10f %16.10f\n",
                        xk, path.yp.host(1, path.kount), std::sin(xk),
                        path.yp.host(2, path.kount), std::cos(xk));
        }

        double err1 = std::fabs(ystart.host(1));
        double err2 = std::fabs(ystart.host(2) - 1.0);
        double max_err = std::max(err1, err2);

        std::printf("\nFinal values at 2*pi:\n");
        std::printf("  y1 = %16.10f  (exact 0)\n", ystart.host(1));
        std::printf("  y2 = %16.10f  (exact 1)\n", ystart.host(2));
        std::printf("  Max error: %.2e\n", max_err);
        std::printf("  Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
