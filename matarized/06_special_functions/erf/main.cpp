#include <cstdio>
#include <cmath>
#include <matar.h>
#include "erf.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 7;
        double test_x[NTEST] = {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0};

        DCArrayKokkos<double> x_d(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x_d.host(i) = test_x[i];
        x_d.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = nr_erf(x_d(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Error function erf(x) via nr_erf:\n");
        std::printf("%12s %18s %18s %14s\n", "x", "nr_erf(x)", "std::erf(x)", "Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double expected = std::erf(test_x[i]);
            double err = fabs(result.host(i) - expected);
            if (err > max_err) max_err = err;
            std::printf("%12.4f %18.10f %18.10f %14.2e\n",
                        test_x[i], result.host(i), expected, err);
        }

        double erf1 = result.host(5);
        double expected_erf1 = 0.8427007929;
        double err1 = fabs(erf1 - expected_erf1);
        std::printf("\nnr_erf(1.0) = %.10f  expected = %.10f  err = %.2e\n",
                    erf1, expected_erf1, err1);
        std::printf("Max error vs std::erf: %.2e\n", max_err);
        std::printf("Test %s\n", err1 < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
