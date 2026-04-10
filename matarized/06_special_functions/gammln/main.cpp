#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gammln.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 10;
        double test_x[NTEST] = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

        DCArrayKokkos<double> x(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x.host(i) = test_x[i];
        x.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = gammln(x(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Log of gamma function (MATAR parallel evaluation):\n");
        std::printf("%12s %18s %18s %14s\n", "X", "GAMMLN(X)", "lgamma(X)", "Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::lgamma(test_x[i]);
            double err = std::fabs(result.host(i) - expected);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n",
                        test_x[i], result.host(i), expected, err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
