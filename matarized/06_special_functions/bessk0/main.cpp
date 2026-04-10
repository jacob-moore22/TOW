#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bessk0.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 6;
        double test_x[NTEST] = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};

        DCArrayKokkos<double> x(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++)
            x.host(i) = test_x[i];
        x.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = bessk0(x(i));
        });
        MATAR_FENCE();
        result.update_host();

        constexpr double ref_K0_1 = 0.4210244382407084;

        std::printf("Modified Bessel K0(x) — MATAR parallel evaluation:\n");
        std::printf("%12s %20s\n", "x", "K0(x)");
        for (int i = 0; i < NTEST; i++)
            std::printf("%12.4f %20.12f\n", test_x[i], result.host(i));

        double err = std::fabs(result.host(2) - ref_K0_1);
        std::printf("\nK0(1.0) = %.16f  (ref %.16f)\n", result.host(2), ref_K0_1);
        std::printf("Error: %.2e\n", err);
        std::printf("Test %s\n", err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
