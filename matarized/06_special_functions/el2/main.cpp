#include <cstdio>
#include <cmath>
#include <matar.h>
#include "el2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // When kc = 1 (k = 0), el2(x, 1, 1, 1) = atan(x)
        constexpr int NTEST = 12;
        double test_x[NTEST] = {
             0.1,  0.5,  1.0,  2.0,   5.0,  10.0,
            -0.1, -0.5, -1.0, -2.0,  -5.0, -10.0
        };

        DCArrayKokkos<double> x_arr(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x_arr.host(i) = test_x[i];
        x_arr.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = el2(x_arr(i), 1.0, 1.0, 1.0);
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("General Elliptic Integral EL2 (MATAR parallel evaluation):\n\n");
        std::printf("Test: el2(x, 1, 1, 1) = atan(x)\n");
        std::printf("%12s %18s %18s %14s\n",
                    "x", "el2(x,1,1,1)", "atan(x)", "Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::atan(test_x[i]);
            double err = std::fabs(result.host(i) - expected);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n",
                        test_x[i], result.host(i), expected, err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        bool pass = (max_err < 1e-6);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
