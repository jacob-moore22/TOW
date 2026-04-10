#include <cstdio>
#include <cmath>
#include <complex>
#include <matar.h>
#include "laguer.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== LAGUER: Laguerre's method for polynomial roots ===\n\n");

        // p(x) = x^2 - 2 = -2 + 0*x + 1*x^2  → roots at +/-sqrt(2)
        {
            std::complex<double> a[3] = {{-2.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}};
            std::complex<double> x(1.0, 0.0);
            laguer(a, 2, x, 1.0e-12, true);
            std::printf("x^2 - 2:  root = (%.15f, %.15f)\n", x.real(), x.imag());
            std::printf("  expected sqrt(2) = %.15f, err = %.2e\n",
                        std::sqrt(2.0), std::fabs(x.real() - std::sqrt(2.0)));
        }

        // p(x) = x^3 - x - 1  → one real root near 1.3247
        {
            std::complex<double> a[4] = {{-1.0, 0.0}, {-1.0, 0.0},
                                          {0.0, 0.0},  {1.0, 0.0}};
            std::complex<double> x(1.5, 0.0);
            laguer(a, 3, x, 1.0e-12, true);
            std::printf("\nx^3-x-1:  root = (%.15f, %.15f)\n", x.real(), x.imag());
            double val = x.real() * x.real() * x.real() - x.real() - 1.0;
            std::printf("  f(root) = %.2e\n", std::fabs(val));
        }

        // p(x) = x^2 + 1  → roots at +/- i
        {
            std::complex<double> a[3] = {{1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}};
            std::complex<double> x(0.5, 0.5);
            laguer(a, 2, x, 1.0e-12, true);
            std::printf("\nx^2 + 1:  root = (%.15f, %.15f)\n", x.real(), x.imag());
            std::printf("  expected (0, +/-1), err_real = %.2e, err_imag = %.2e\n",
                        std::fabs(x.real()), std::fabs(std::fabs(x.imag()) - 1.0));
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
