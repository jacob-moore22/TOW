#include <cstdio>
#include <cmath>
#include <matar.h>

using namespace mtr;

#include "../../02_linear_algebra/ludcmp/ludcmp.hpp"
#include "../../02_linear_algebra/lubksb/lubksb.hpp"
#include "../../16_ode_integration/rk4/rk4.hpp"
#include "../../16_ode_integration/rkqc/rkqc.hpp"
#include "../../16_ode_integration/odeint/odeint.hpp"
#include "shootf.hpp"

// Test: y'' + y = 0, y(0)=0, y(pi)=0
// Two-sided shooting to xf = pi/2.
// Left: y(0)=0, free v1(1)=y'(0)
// Right: y(pi)=0, free v2(1)=y'(pi)
// Exact: y=sin(x), y'(0)=1, y'(pi)=-1

void derivs_test(double /*x*/, DFMatrixKokkos<double>& y,
                 DFMatrixKokkos<double>& dydx)
{
    dydx.host(1) =  y.host(2);
    dydx.host(2) = -y.host(1);
}

void load1_test(double /*x1*/, DFMatrixKokkos<double>& v1,
                DFMatrixKokkos<double>& y)
{
    y.host(1) = 0.0;
    y.host(2) = v1.host(1);
}

void load2_test(double /*x2*/, DFMatrixKokkos<double>& v2,
                DFMatrixKokkos<double>& y)
{
    y.host(1) = 0.0;
    y.host(2) = v2.host(1);
}

void score_test(double /*xf*/, DFMatrixKokkos<double>& y,
                DFMatrixKokkos<double>& f)
{
    f.host(1) = y.host(1);
    f.host(2) = y.host(2);
}

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NVAR = 2;
        constexpr int N1   = 1;
        constexpr int N2   = 1;

        DFMatrixKokkos<double> v1(N2), v2(N1);
        DFMatrixKokkos<double> delv1(N2), delv2(N1);
        DFMatrixKokkos<double> f(NVAR), dv1(N2), dv2(N1);

        v1.host(1) = 0.5;   delv1.host(1) = 0.01;
        v2.host(1) = -0.5;  delv2.host(1) = 0.01;

        double x1 = 0.0;
        double x2 = M_PI;
        double xf = M_PI / 2.0;

        std::printf("shootf test: y'' + y = 0, y(0)=0, y(pi)=0, fit at pi/2\n");
        std::printf("Initial: v1(1)=y'(0)=%.4f, v2(1)=y'(pi)=%.4f\n",
                    v1.host(1), v2.host(1));

        auto stepper = [](DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                          int n, double& x, double htry, double eps,
                          DFMatrixKokkos<double>& yscal,
                          double& hdid, double& hnext, auto derivs) {
            rkqc(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs);
        };

        for (int iter = 0; iter < 10; iter++) {
            shootf(NVAR, v1, v2, delv1, delv2, N1, N2,
                   x1, x2, xf, 1.0e-6, 0.1, 1.0e-10,
                   f, dv1, dv2,
                   load1_test, load2_test, score_test, derivs_test, stepper);
            std::printf("  iter %d: y'(0)=%10.6f  y'(pi)=%10.6f\n",
                        iter + 1, v1.host(1), v2.host(1));
            if (std::fabs(v1.host(1) - 1.0) < 1e-8 &&
                std::fabs(v2.host(1) + 1.0) < 1e-8)
                break;
        }

        double err1 = std::fabs(v1.host(1) - 1.0);
        double err2 = std::fabs(v2.host(1) + 1.0);
        std::printf("Final: y'(0)=%.10f (exact 1.0, err=%.2e)\n", v1.host(1), err1);
        std::printf("       y'(pi)=%.10f (exact -1.0, err=%.2e)\n", v2.host(1), err2);
        std::printf("Test %s\n", (err1 < 1e-4 && err2 < 1e-4) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
