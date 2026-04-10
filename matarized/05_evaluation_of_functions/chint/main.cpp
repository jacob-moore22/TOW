#include <cstdio>
#include <cmath>
#include <matar.h>
#include "../chebft/chebft.hpp"
#include "../chebev/chebev.hpp"
#include "chint.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;
        double a = -1.0, b = 1.0;

        DFMatrixKokkos<double> c(N);
        DFMatrixKokkos<double> cint(N);

        chebft(a, b, c, N, [](double x) { return std::exp(x); });
        chint(a, b, c, cint, N);

        std::printf("Chebyshev integral of exp(x) on [-1, 1]:\n");
        std::printf("(integral from -1 to x of exp(t) dt = exp(x) - exp(-1))\n\n");
        std::printf("%12s %18s %18s %14s\n", "X", "Exact", "chint+chebev", "Error");

        double max_err = 0.0;
        double base = std::exp(-1.0);
        for (int i = -10; i <= 10; i += 2) {
            double x     = i * 0.1;
            double integ = chebev(a, b, cint, N, x);
            double exact = std::exp(x) - base;
            double err   = std::fabs(integ - exact);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n", x, exact, integ, err);
        }
        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-8 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
