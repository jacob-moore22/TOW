#include <cstdio>
#include <cmath>
#include <matar.h>
#include "erfcc.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 13;
        double test_x[NTEST] = {
            -3.0, -2.0, -1.0, -0.5, 0.0, 0.25, 0.5,
             1.0,  1.5,  2.0,  2.5,  3.0,  4.0
        };

        DCArrayKokkos<double> x_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x_vals.host(i) = test_x[i];
        x_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = erfcc(x_vals(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Complementary error function (MATAR parallel evaluation):\n");
        std::printf("%10s %18s %18s %14s\n", "X", "ERFCC(X)", "erfc(X)", "Abs Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::erfc(test_x[i]);
            double err = std::fabs(result.host(i) - expected);
            if (err > max_err) max_err = err;
            std::printf("%10.4f %18.10f %18.10f %14.2e\n",
                        test_x[i], result.host(i), expected, err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
