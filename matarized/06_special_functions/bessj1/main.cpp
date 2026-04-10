#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bessj1.hpp"

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
            result(i) = bessj1(x(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Bessel J1(x) — MATAR parallel evaluation:\n");
        std::printf("%12s %18s %18s %14s\n", "X", "BESSJ1(X)", "std::ref", "Error");

        double ref[NTEST] = {
            0.2422684576748739, 0.4400505857449335, 0.5767248077568734,
            0.3390589585259365, -0.3275791375914652, 0.2346363468539146,
            0.0434727461688614, 0.0668331429926655
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
