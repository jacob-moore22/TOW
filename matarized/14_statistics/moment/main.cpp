#include <cstdio>
#include <cmath>
#include <matar.h>
#include "moment.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 1000;
        DFMatrixKokkos<double> data(N);

        // Fill with linearly spaced values 1..N so statistics are known
        for (int i = 1; i <= N; i++)
            data.host(i) = static_cast<double>(i);

        double ave, adev, sdev, var, skew, curt;
        moment(data, N, ave, adev, sdev, var, skew, curt);

        // Expected for uniform 1..N:
        //   mean     = (N+1)/2
        //   variance = (N^2-1)/12
        double exp_ave = (N + 1) / 2.0;
        double exp_var = (static_cast<double>(N) * N - 1.0) / 12.0;

        std::printf("Moment statistics for data = 1..%d:\n", N);
        std::printf("  Mean     = %14.6f  (expected %14.6f)\n", ave, exp_ave);
        std::printf("  Avg dev  = %14.6f\n", adev);
        std::printf("  Std dev  = %14.6f  (expected %14.6f)\n", sdev, sqrt(exp_var));
        std::printf("  Variance = %14.6f  (expected %14.6f)\n", var, exp_var);
        std::printf("  Skewness = %14.6f  (expected %14.6f)\n", skew, 0.0);
        std::printf("  Kurtosis = %14.6f  (expected %14.6f)\n", curt, -1.2);

        double tol = 1e-6;
        bool pass = fabs(ave - exp_ave) < tol &&
                    fabs(var - exp_var) / exp_var < tol;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
