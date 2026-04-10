#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bessj0.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 8;
        double test_x[NTEST] = {0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 20.0};

        DCArrayKokkos<double> x(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x.host(i) = test_x[i];
        x.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = bessj0(x(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Bessel J0(x) — MATAR parallel evaluation:\n");
        std::printf("%12s %18s %18s %14s\n", "X", "BESSJ0(X)", "std::ref", "Error");

        double ref[NTEST] = {
            0.9384698072408129, 0.7651976865579666, 0.2238907791412357,
           -0.2600519549019334, -0.1775967713143383, 0.1716508071375539,
           -0.2459357644513483, 0.1670246643401472
        };

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double err = std::fabs(result.host(i) - ref[i]);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n",
                        test_x[i], result.host(i), ref[i], err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
