#include <cstdio>
#include <cmath>
#include <matar.h>
#include "beta.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 8;
        double test_z[NTEST] = {0.5, 1.0, 1.0, 2.0, 2.0, 3.0,  5.0,  10.0};
        double test_w[NTEST] = {0.5, 1.0, 2.0, 2.0, 3.0, 4.0, 10.0,  10.0};

        DCArrayKokkos<double> z_vals(NTEST);
        DCArrayKokkos<double> w_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++) {
            z_vals.host(i) = test_z[i];
            w_vals.host(i) = test_w[i];
        }
        z_vals.update_device();
        w_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = beta(z_vals(i), w_vals(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Beta function (MATAR parallel evaluation):\n");
        std::printf("%8s %8s %18s %18s %14s\n",
                    "Z", "W", "BETA(Z,W)", "Reference", "Rel Error");

        double max_rel_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::exp(
                std::lgamma(test_z[i]) + std::lgamma(test_w[i])
                - std::lgamma(test_z[i] + test_w[i]));
            double rel_err = (expected != 0.0)
                           ? std::fabs(result.host(i) - expected) / expected
                           : std::fabs(result.host(i));
            if (rel_err > max_rel_err) max_rel_err = rel_err;
            std::printf("%8.2f %8.2f %18.10f %18.10f %14.2e\n",
                        test_z[i], test_w[i], result.host(i), expected, rel_err);
        }

        std::printf("\nMax relative error: %.2e\n", max_rel_err);
        std::printf("Test %s\n", max_rel_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
