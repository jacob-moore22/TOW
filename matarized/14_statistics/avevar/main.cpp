#include <cstdio>
#include <cmath>
#include <matar.h>
#include "avevar.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 500;
        DFMatrixKokkos<double> data(N);

        for (int i = 1; i <= N; i++)
            data.host(i) = static_cast<double>(i);

        double ave, var;
        avevar(data, N, ave, var);

        double exp_ave = (N + 1) / 2.0;
        double exp_var = (static_cast<double>(N) * N - 1.0) / 12.0;

        std::printf("Avevar for data = 1..%d:\n", N);
        std::printf("  Mean     = %14.6f  (expected %14.6f)\n", ave, exp_ave);
        std::printf("  Variance = %14.6f  (expected %14.6f)\n", var, exp_var);

        bool pass = fabs(ave - exp_ave) < 1e-10 &&
                    fabs(var - exp_var) / exp_var < 1e-10;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
