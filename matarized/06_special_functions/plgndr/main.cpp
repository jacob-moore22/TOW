#include <cstdio>
#include <cmath>
#include <matar.h>
#include "plgndr.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Analytic reference values at x = 0.5:
        //   P_0^0 = 1
        //   P_1^0 = x = 0.5
        //   P_1^1 = -sqrt(1-x^2) = -sqrt(0.75)
        //   P_2^0 = (3x^2-1)/2 = -0.125
        //   P_2^1 = -3x*sqrt(1-x^2)
        //   P_2^2 = 3(1-x^2) = 2.25
        //   P_3^0 = (5x^3-3x)/2 = -0.4375
        //   P_3^1 = -1.5*(5x^2-1)*sqrt(1-x^2)
        //   P_4^0 = (35x^4-30x^2+3)/8
        //   P_5^0 = (63x^5-70x^3+15x)/8

        constexpr double xv = 0.5;
        double sx = std::sqrt(1.0 - xv * xv);

        constexpr int NTEST = 10;
        int test_l[NTEST] = {0, 1, 1, 2, 2, 2, 3, 3, 4, 5};
        int test_m[NTEST] = {0, 0, 1, 0, 1, 2, 0, 1, 0, 0};
        double ref[NTEST] = {
            1.0,
            xv,
            -sx,
            (3.0 * xv * xv - 1.0) / 2.0,
            -3.0 * xv * sx,
            3.0 * (1.0 - xv * xv),
            (5.0 * xv * xv * xv - 3.0 * xv) / 2.0,
            -1.5 * (5.0 * xv * xv - 1.0) * sx,
            (35.0 * xv * xv * xv * xv - 30.0 * xv * xv + 3.0) / 8.0,
            (63.0 * xv * xv * xv * xv * xv - 70.0 * xv * xv * xv + 15.0 * xv) / 8.0
        };

        DCArrayKokkos<int>    l_vals(NTEST);
        DCArrayKokkos<int>    m_vals(NTEST);
        DCArrayKokkos<double> result(NTEST);

        for (int i = 0; i < NTEST; i++) {
            l_vals.host(i) = test_l[i];
            m_vals.host(i) = test_m[i];
        }
        l_vals.update_device();
        m_vals.update_device();

        FOR_ALL(i, 0, NTEST, {
            result(i) = plgndr(l_vals(i), m_vals(i), xv);
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Associated Legendre polynomials at x = %.2f "
                    "(MATAR parallel evaluation):\n", xv);
        std::printf("%4s %4s %18s %18s %14s\n",
                    "L", "M", "PLGNDR", "Reference", "Abs Error");

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double err = std::fabs(result.host(i) - ref[i]);
            if (err > max_err) max_err = err;
            std::printf("%4d %4d %18.10f %18.10f %14.2e\n",
                        test_l[i], test_m[i], result.host(i), ref[i], err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-14 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
