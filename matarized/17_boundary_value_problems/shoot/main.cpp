#include <cstdio>
#include <cmath>
#include <matar.h>

using namespace mtr;

#include "../../02_linear_algebra/ludcmp/ludcmp.hpp"
#include "../../02_linear_algebra/lubksb/lubksb.hpp"
#include "../../16_ode_integration/rk4/rk4.hpp"
#include "../../16_ode_integration/rkqc/rkqc.hpp"
#include "../../16_ode_integration/odeint/odeint.hpp"
#include "shoot.hpp"

// Test problem: y'' + y = 0, y(0)=0, y(pi/2)=1
// Exact: y = sin(x), so y'(0) = 1.
// NVAR=2: y[1]=y, y[2]=y'
// Free parameter: v[1] = y'(0)

void derivs_test(double /*x*/, DFMatrixKokkos<double>& y,
                 DFMatrixKokkos<double>& dydx)
{
    dydx.host(1) =  y.host(2);
    dydx.host(2) = -y.host(1);
}

void load_test(double /*x1*/, DFMatrixKokkos<double>& v,
               DFMatrixKokkos<double>& y)
{
    y.host(1) = 0.0;
    y.host(2) = v.host(1);
}

void score_test(double /*x2*/, DFMatrixKokkos<double>& y,
                DFMatrixKokkos<double>& f)
{
    f.host(1) = y.host(1) - 1.0;
}

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NVAR = 2;
        constexpr int N2   = 1;

        DFMatrixKokkos<double> v(N2);
        DFMatrixKokkos<double> delv(N2);
        DFMatrixKokkos<double> f(N2);
        DFMatrixKokkos<double> dv_arr(N2);

        v.host(1)    = 0.5;
        delv.host(1) = 0.01;

        double x1 = 0.0;
        double x2 = M_PI / 2.0;

        std::printf("Shooting method test: y'' + y = 0, y(0)=0, y(pi/2)=1\n");
        std::printf("Initial guess: y'(0) = %.4f\n", v.host(1));

        auto stepper = [](DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                          int n, double& x, double htry, double eps,
                          DFMatrixKokkos<double>& yscal,
                          double& hdid, double& hnext, auto derivs) {
            rkqc(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs);
        };

        for (int iter = 0; iter < 10; iter++) {
            shoot(NVAR, v, delv, N2, x1, x2, 1.0e-6, 0.1, 1.0e-10,
                  f, dv_arr, load_test, score_test, derivs_test, stepper);
            std::printf("  iter %d: v[1]=y'(0) = %.10f\n", iter + 1, v.host(1));
            if (std::fabs(v.host(1) - 1.0) < 1.0e-8) break;
        }

        double error = std::fabs(v.host(1) - 1.0);
        std::printf("Final y'(0) = %.10f  (exact=1.0, error=%.2e)\n", v.host(1), error);
        std::printf("Test %s\n", error < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
