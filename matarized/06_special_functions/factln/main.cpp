#include <cstdio>
#include <cmath>
#include <matar.h>
#include "factln.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 14;
        int test_n[NTEST] = {0, 1, 2, 3, 4, 5, 10, 15, 20, 50, 100, 150, 200, 500};

        DCArrayKokkos<int>    n_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            n_vals.host(i) = test_n[i];
        n_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = factln(n_vals(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Log factorial function (MATAR parallel evaluation):\n");
        std::printf("%6s %18s %18s %14s\n", "N", "FACTLN(N)", "lgamma(N+1)", "Abs Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::lgamma(test_n[i] + 1.0);
            double err = std::fabs(result.host(i) - expected);
            if (err > max_err) max_err = err;
            std::printf("%6d %18.10f %18.10f %14.2e\n",
                        test_n[i], result.host(i), expected, err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
