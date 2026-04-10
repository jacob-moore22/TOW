#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ratint.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 5;
        DFMatrixKokkos<double> xa(N);
        DFMatrixKokkos<double> ya(N);

        for (int i = 1; i <= N; i++) {
            double xi = static_cast<double>(i - 1);
            xa.host(i) = xi;
            ya.host(i) = 1.0 / (1.0 + xi * xi);  // f(x) = 1/(1+x^2)
        }

        double x_test[] = {0.5, 1.5, 2.5, 3.5};
        constexpr int NTEST = 4;

        std::printf("Rational interpolation of f(x) = 1/(1+x^2) using %d points\n", N);
        std::printf("%10s %16s %16s %16s %14s\n",
                    "x", "ratint(x)", "exact", "dy (err est)", "actual err");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double y, dy;
            ratint(xa, ya, N, x_test[t], y, dy);
            double exact = 1.0 / (1.0 + x_test[t] * x_test[t]);
            double err   = std::fabs(y - exact);
            if (err > max_err) max_err = err;
            std::printf("%10.4f %16.10f %16.10f %16.2e %14.2e\n",
                        x_test[t], y, exact, dy, err);
        }

        std::printf("\nMax actual error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-2 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
