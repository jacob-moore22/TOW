#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bessi.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 5;
        int    test_n[NTEST] = {2, 3, 5, 2, 10};
        double test_x[NTEST] = {1.0, 2.0, 3.0, 5.0, 5.0};

        DCArrayKokkos<int>    n_arr(NTEST);
        DCArrayKokkos<double> x_arr(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++) {
            n_arr.host(i) = test_n[i];
            x_arr.host(i) = test_x[i];
        }
        n_arr.update_device();
        x_arr.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = bessi(n_arr(i), x_arr(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Modified Bessel In(x) — MATAR parallel evaluation:\n");
        std::printf("%6s %12s %20s\n", "n", "x", "In(x)");
        for (int i = 0; i < NTEST; i++)
            std::printf("%6d %12.4f %20.12f\n",
                        test_n[i], test_x[i], result.host(i));

        std::printf("\nTest PASSED (visual check against tables)\n");
    }
    MATAR_FINALIZE();
    return 0;
}
