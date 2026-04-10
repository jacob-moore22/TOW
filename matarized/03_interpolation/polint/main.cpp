#include <cstdio>
#include <cmath>
#include <matar.h>
#include "polint.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 5;
        DFMatrixKokkos<double> xa(N);
        DFMatrixKokkos<double> ya(N);

        for (int i = 1; i <= N; i++) {
            xa.host(i) = static_cast<double>(i);
            ya.host(i) = xa.host(i) * xa.host(i) * xa.host(i);  // y = x^3
        }

        double x_test[] = {1.5, 2.5, 3.5, 4.5};
        constexpr int NTEST = 4;

        std::printf("Polynomial interpolation of f(x) = x^3 using %d points\n", N);
        std::printf("%10s %16s %16s %16s %14s\n",
                    "x", "polint(x)", "exact", "dy (err est)", "actual err");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double y, dy;
            polint(xa, ya, N, x_test[t], y, dy);
            double exact = x_test[t] * x_test[t] * x_test[t];
            double err   = std::fabs(y - exact);
            if (err > max_err) max_err = err;
            std::printf("%10.4f %16.10f %16.10f %16.2e %14.2e\n",
                        x_test[t], y, exact, dy, err);
        }

        std::printf("\nMax actual error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
