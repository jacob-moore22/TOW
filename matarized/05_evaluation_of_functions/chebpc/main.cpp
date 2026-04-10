#include <cstdio>
#include <cmath>
#include <matar.h>
#include "../chebft/chebft.hpp"
#include "chebpc.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;
        double a = -1.0, b = 1.0;

        DFMatrixKokkos<double> c(N);
        DFMatrixKokkos<double> d(N);

        chebft(a, b, c, N, [](double x) { return std::exp(x); });
        chebpc(c, d, N);

        std::printf("Chebyshev-to-polynomial for exp(x) on [-1, 1]:\n\n");
        std::printf("%12s %18s %18s %14s\n", "X", "exp(X)", "Polynomial", "Error");

        double max_err = 0.0;
        for (int i = -10; i <= 10; i += 2) {
            double x = i * 0.1;
            double y = (2.0 * x - a - b) / (b - a);

            double poly = d(N);
            for (int j = N - 1; j >= 1; j--)
                poly = poly * y + d(j);

            double exact = std::exp(x);
            double err   = std::fabs(poly - exact);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n", x, exact, poly, err);
        }
        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-8 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
