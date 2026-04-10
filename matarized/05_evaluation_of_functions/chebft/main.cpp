#include <cstdio>
#include <cmath>
#include <matar.h>
#include "chebft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;
        double a = -1.0, b = 1.0;

        DFMatrixKokkos<double> c(N);
        chebft(a, b, c, N, [](double x) { return std::exp(x); });

        std::printf("Chebyshev coefficients for exp(x) on [-1, 1]:\n");
        for (int j = 1; j <= N; j++)
            std::printf("  c(%2d) = %20.12e\n", j, c(j));

        std::printf("\n%12s %18s %18s %14s\n", "X", "exp(X)", "Cheb fit", "Error");
        double max_err = 0.0;
        for (int i = -10; i <= 10; i += 2) {
            double x = i * 0.1;
            double y = (2.0 * x - a - b) / (b - a);

            double d = 0.0, dd = 0.0;
            for (int j = N; j >= 2; j--) {
                double sv = d;
                d = 2.0 * y * d - dd + c(j);
                dd = sv;
            }
            double fit = y * d - dd + 0.5 * c(1);
            double exact = std::exp(x);
            double err = std::fabs(fit - exact);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n", x, exact, fit, err);
        }
        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
