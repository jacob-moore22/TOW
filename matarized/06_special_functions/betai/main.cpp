#include <cstdio>
#include <cmath>
#include <matar.h>
#include "betai.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 5;
        double test_a[NTEST] = {1.0, 1.0, 2.0, 0.5, 5.0};
        double test_b[NTEST] = {1.0, 1.0, 3.0, 0.5, 5.0};
        double test_x[NTEST] = {0.5, 0.25, 0.4, 0.5, 0.5};

        DCArrayKokkos<double> a_d(NTEST);
        DCArrayKokkos<double> b_d(NTEST);
        DCArrayKokkos<double> x_d(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++) {
            a_d.host(i) = test_a[i];
            b_d.host(i) = test_b[i];
            x_d.host(i) = test_x[i];
        }
        a_d.update_device();
        b_d.update_device();
        x_d.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = betai(a_d(i), b_d(i), x_d(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Incomplete beta function I_x(a,b) via betai:\n");
        std::printf("%8s %8s %8s %18s\n", "a", "b", "x", "I_x(a,b)");
        for (int i = 0; i < NTEST; i++) {
            std::printf("%8.2f %8.2f %8.2f %18.10f\n",
                        test_a[i], test_b[i], test_x[i], result.host(i));
        }

        double expected = 0.5;
        double err = fabs(result.host(0) - expected);
        std::printf("\nbetai(1,1,0.5) = %.10f  expected = %.10f  err = %.2e\n",
                    result.host(0), expected, err);
        std::printf("Test %s\n", err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
