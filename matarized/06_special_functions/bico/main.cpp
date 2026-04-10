#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bico.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 10;
        int test_n[NTEST] = { 1,  4,  5, 10, 10, 20, 20, 30, 50, 100};
        int test_k[NTEST] = { 0,  2,  5,  3,  5, 10,  7, 15, 25,  50};
        double ref[NTEST] = {
            1.0, 6.0, 1.0, 120.0, 252.0,
            184756.0, 77520.0, 155117520.0, 1.26410606437752e14,
            1.00891344545564e29
        };

        DCArrayKokkos<int>    n_vals(NTEST);
        DCArrayKokkos<int>    k_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++) {
            n_vals.host(i) = test_n[i];
            k_vals.host(i) = test_k[i];
        }
        n_vals.update_device();
        k_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = bico(n_vals(i), k_vals(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Binomial coefficient (MATAR parallel evaluation):\n");
        std::printf("%6s %6s %22s %22s %14s\n",
                    "N", "K", "BICO(N,K)", "Reference", "Rel Error");

        double max_rel_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double rel_err = (ref[i] != 0.0)
                           ? std::fabs(result.host(i) - ref[i]) / ref[i]
                           : std::fabs(result.host(i));
            if (rel_err > max_rel_err) max_rel_err = rel_err;
            std::printf("%6d %6d %22.6f %22.6f %14.2e\n",
                        test_n[i], test_k[i], result.host(i), ref[i], rel_err);
        }

        std::printf("\nMax relative error: %.2e\n", max_rel_err);
        std::printf("Test %s\n", max_rel_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
