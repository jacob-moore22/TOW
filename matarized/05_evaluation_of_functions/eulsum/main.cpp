#include <cstdio>
#include <cmath>
#include <matar.h>
#include "eulsum.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NVAL = 40;

        std::printf("Euler summation of ln(1+x) = x - x^2/2 + x^3/3 - ...\n\n");
        std::printf("%12s %18s %18s %14s\n", "X", "ln(1+X)", "Euler sum", "Error");

        DFMatrixKokkos<double> wksp(NVAL);

        double max_err = 0.0;
        for (int i = -8; i <= 8; i++) {
            double x = i / 10.0;

            for (int j = 1; j <= NVAL; j++) wksp(j) = 0.0;

            double sum    = 0.0;
            double xpower = -1.0;
            int    nterm  = 0;

            for (int j = 1; j <= NVAL; j++) {
                xpower = -x * xpower;
                double term = xpower / j;
                eulsum(sum, term, j, wksp, nterm);
            }

            double exact = std::log(1.0 + x);
            double err   = std::fabs(sum - exact);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n", x, exact, sum, err);
        }
        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
