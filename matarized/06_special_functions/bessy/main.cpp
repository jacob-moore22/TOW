#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bessy.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 6;
        int    test_n[NTEST] = {2,    3,    5,    2,     3,     10};
        double test_x[NTEST] = {1.0,  2.0,  5.0,  10.0,  10.0,  10.0};

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
            result(i) = bessy(n_arr(i), x_arr(i));
        });
        MATAR_FENCE();
        result.update_host();

        std::printf("Bessel Yn(x) — MATAR parallel evaluation:\n");
        std::printf("%6s %12s %18s %18s %14s\n", "N", "X", "BESSY(N,X)", "std::ref", "Error");

        double ref[NTEST] = {
           -1.6506826068162546,   // Y2(1)
           -1.1277837768404277,   // Y3(2)
            0.4536948224911019,   // Y5(5)
           -0.0058680824422086,   // Y2(10)
            0.2511895880110441,   // Y3(10)
           -0.3598141521834027    // Y10(10)
        };

        double max_err = 0.0;
        for (int i = 0; i < NTEST; i++) {
            double err = std::fabs(result.host(i) - ref[i]);
            if (err > max_err) max_err = err;
            std::printf("%6d %12.4f %18.10f %18.10f %14.2e\n",
                        test_n[i], test_x[i], result.host(i), ref[i], err);
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-4 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
