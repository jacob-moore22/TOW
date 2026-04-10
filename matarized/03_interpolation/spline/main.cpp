#include <cstdio>
#include <cmath>
#include <matar.h>
#include "spline.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 10;
        DFMatrixKokkos<double> x(N);
        DFMatrixKokkos<double> y(N);
        DFMatrixKokkos<double> y2(N);

        for (int i = 1; i <= N; i++) {
            x.host(i) = (i - 1) * M_PI / (N - 1);
            y.host(i) = std::sin(x.host(i));
        }

        spline(x, y, N, std::cos(0.0), std::cos(M_PI), y2);

        std::printf("Cubic spline second derivatives for sin(x) on [0, pi]\n");
        std::printf("%6s %14s %14s %14s\n", "i", "x(i)", "y(i)", "y2(i)");
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %14.8f %14.8f %14.8f\n",
                        i, x.host(i), y.host(i), y2.host(i));
        }

        // For sin(x), y'' = -sin(x); verify at interior points
        std::printf("\nVerification: y2 vs -sin(x)\n");
        std::printf("%6s %14s %14s %14s\n", "i", "y2(i)", "-sin(x)", "error");
        double max_err = 0.0;
        for (int i = 2; i <= N - 1; i++) {
            double exact = -std::sin(x.host(i));
            double err   = std::fabs(y2.host(i) - exact);
            if (err > max_err) max_err = err;
            std::printf("%6d %14.8f %14.8f %14.2e\n",
                        i, y2.host(i), exact, err);
        }
        std::printf("\nMax error at interior points: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 0.05 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
