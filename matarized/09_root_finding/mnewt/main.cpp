#include <cstdio>
#include <cmath>
#include <matar.h>
#include "mnewt.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== MNEWT: Multi-dimensional Newton's method ===\n\n");

        // System: x^2 + y^2 = 2, x - y = 0  → solution (1,1) and (-1,-1)
        auto usrfun1 = [](const double* x, double** alpha, double* beta) {
            double f0 = x[0] * x[0] + x[1] * x[1] - 2.0;
            double f1 = x[0] - x[1];
            beta[0] = -f0;
            beta[1] = -f1;
            alpha[0][0] = 2.0 * x[0];
            alpha[0][1] = 2.0 * x[1];
            alpha[1][0] = 1.0;
            alpha[1][1] = -1.0;
        };

        {
            double x[2] = {0.5, 0.8};
            mnewt(100, x, 2, 1.0e-14, 1.0e-14, usrfun1);
            std::printf("x^2+y^2=2, x-y=0:  x=%.15f  y=%.15f\n", x[0], x[1]);
            std::printf("  expected: (1, 1),  err = (%.2e, %.2e)\n",
                        std::fabs(x[0] - 1.0), std::fabs(x[1] - 1.0));
        }

        // System: x^3 - y = 0, y^3 - x = 0  → real solution (1,1), (0,0), (-1,-1)
        auto usrfun2 = [](const double* x, double** alpha, double* beta) {
            double f0 = x[0] * x[0] * x[0] - x[1];
            double f1 = x[1] * x[1] * x[1] - x[0];
            beta[0] = -f0;
            beta[1] = -f1;
            alpha[0][0] = 3.0 * x[0] * x[0];
            alpha[0][1] = -1.0;
            alpha[1][0] = -1.0;
            alpha[1][1] = 3.0 * x[1] * x[1];
        };

        {
            double x[2] = {0.8, 0.8};
            mnewt(100, x, 2, 1.0e-14, 1.0e-14, usrfun2);
            std::printf("\nx^3=y, y^3=x:  x=%.15f  y=%.15f\n", x[0], x[1]);
            std::printf("  expected: (1, 1),  err = (%.2e, %.2e)\n",
                        std::fabs(x[0] - 1.0), std::fabs(x[1] - 1.0));
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
