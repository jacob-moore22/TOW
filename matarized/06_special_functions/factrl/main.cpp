#include <cstdio>
#include <cmath>
#include <matar.h>
#include "factrl.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 12;
        int test_n[NTEST] = {0, 1, 2, 3, 4, 5, 6, 10, 15, 20, 25, 32};

        DCArrayKokkos<int>    n_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            n_vals.host(i) = test_n[i];
        n_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = factrl(n_vals(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Factorial function (MATAR parallel evaluation):\n");
        std::printf("%6s %22s %22s %14s\n", "N", "FACTRL(N)", "tgamma(N+1)", "Rel Error");

        double max_rel_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::tgamma(test_n[i] + 1.0);
            double rel_err  = (expected != 0.0)
                            ? std::fabs(result.host(i) - expected) / expected
                            : std::fabs(result.host(i));
            if (rel_err > max_rel_err) max_rel_err = rel_err;
            std::printf("%6d %22.15e %22.15e %14.2e\n",
                        test_n[i], result.host(i), expected, rel_err);
        }

        // Also test large-n path (n > 32, uses gammln)
        constexpr int NLARGE = 3;
        int large_n[NLARGE] = {50, 100, 170};

        DCArrayKokkos<int>    n_large(NLARGE);
        DCArrayKokkos<double> r_large(NLARGE);

        for (int i = 0; i < NLARGE; i++)
            n_large.host(i) = large_n[i];
        n_large.update_device();

        FOR_ALL(i, 0, NLARGE, {
            r_large(i) = factrl(n_large(i));
        });
        MATAR_FENCE();
        r_large.update_host();

        std::printf("\nLarge-n test (gammln path):\n");
        std::printf("%6s %22s %22s %14s\n", "N", "FACTRL(N)", "tgamma(N+1)", "Rel Error");

        for (int i = 0; i < NLARGE; i++) {
            double expected = std::tgamma(large_n[i] + 1.0);
            double rel_err  = (expected != 0.0 && std::isfinite(expected))
                            ? std::fabs(r_large.host(i) - expected) / expected
                            : 0.0;
            if (rel_err > max_rel_err) max_rel_err = rel_err;
            std::printf("%6d %22.15e %22.15e %14.2e\n",
                        large_n[i], r_large.host(i), expected, rel_err);
        }

        std::printf("\nMax relative error: %.2e\n", max_rel_err);
        std::printf("Test %s\n", max_rel_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
