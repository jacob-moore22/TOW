#include <cstdio>
#include <cmath>
#include <matar.h>
#include "../chebft/chebft.hpp"
#include "../chebpc/chebpc.hpp"
#include "pcshft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;
        double a = 0.0, b = 2.0;

        DFMatrixKokkos<double> c(N);
        DFMatrixKokkos<double> d(N);

        chebft(a, b, c, N, [](double x) { return std::exp(x); });
        chebpc(c, d, N);
        pcshft(a, b, d, N);

        std::printf("Shifted polynomial for exp(x) on [0, 2]:\n\n");
        std::printf("%12s %18s %18s %14s\n", "X", "exp(X)", "Polynomial", "Error");

        double max_err = 0.0;
        for (int i = 0; i <= 20; i += 2) {
            double x = i * 0.1;

            double poly = d(N);
            for (int j = N - 1; j >= 1; j--)
                poly = poly * x + d(j);

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
