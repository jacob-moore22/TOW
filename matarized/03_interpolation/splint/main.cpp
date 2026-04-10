#include <cstdio>
#include <cmath>
#include <matar.h>
#include "spline.hpp"
#include "splint.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 10;
        DFMatrixKokkos<double> xa(N);
        DFMatrixKokkos<double> ya(N);
        DFMatrixKokkos<double> y2a(N);

        for (int i = 1; i <= N; i++) {
            xa.host(i) = (i - 1) * M_PI / (N - 1);
            ya.host(i) = std::sin(xa.host(i));
        }

        spline(xa, ya, N, std::cos(0.0), std::cos(M_PI), y2a);

        constexpr int NTEST = 5;
        double x_test[] = {0.3, 0.7, 1.2, 2.0, 2.8};

        std::printf("Cubic spline interpolation of sin(x)\n");
        std::printf("%10s %16s %16s %14s\n", "x", "splint(x)", "sin(x)", "error");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double y;
            splint(xa, ya, y2a, N, x_test[t], y);
            double exact = std::sin(x_test[t]);
            double err   = std::fabs(y - exact);
            if (err > max_err) max_err = err;
            std::printf("%10.4f %16.10f %16.10f %14.2e\n",
                        x_test[t], y, exact, err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-4 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
