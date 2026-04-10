#include <cstdio>
#include <cmath>
#include <matar.h>
#include "cel.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr double PIO2 = 1.5707963268;

        // -----------------------------------------------------------
        // Test 1: cel(1, 1, 1, 1) = K(0) = pi/2
        // -----------------------------------------------------------
        DCArrayKokkos<double> exact_result(1);
        FOR_ALL(i, 0, 1, {
            exact_result(i) = cel(1.0, 1.0, 1.0, 1.0);
        });
        MATAR_FENCE();
        exact_result.update_host();

        double exact_err = std::fabs(exact_result.host(0) - PIO2);

        std::printf("Complete Elliptic Integral CEL (MATAR parallel evaluation):\n\n");
        std::printf("Test 1: cel(1, 1, 1, 1) = K(0) = pi/2\n");
        std::printf("  Expected: %.10f\n", PIO2);
        std::printf("  Computed: %.10f\n", exact_result.host(0));
        std::printf("  Error:    %.2e\n\n", exact_err);

        // -----------------------------------------------------------
        // Test 2: Legendre relation  E*K' + E'*K - K*K' = pi/2
        //   K(k)  = cel(kc, 1, 1, 1)       with kc = sqrt(1-k^2)
        //   E(k)  = cel(kc, 1, 1, kc^2)
        //   K(k') = cel(k, 1, 1, 1)
        //   E(k') = cel(k, 1, 1, k^2)
        // -----------------------------------------------------------
        constexpr int NLEG = 5;
        double test_k[NLEG] = {0.2, 0.4, 0.5, 0.6, 0.8};

        DCArrayKokkos<double> k_arr(NLEG);
        DCArrayKokkos<double> Kk(NLEG), Ek(NLEG), Kkp(NLEG), Ekp(NLEG);

        for (int i = 0; i < NLEG; i++)
            k_arr.host(i) = test_k[i];
        k_arr.update_device();

        FOR_ALL(i, 0, NLEG, {
            double k  = k_arr(i);
            double kc = sqrt(1.0 - k * k);
            Kk(i)  = cel(kc, 1.0, 1.0, 1.0);
            Ek(i)  = cel(kc, 1.0, 1.0, kc * kc);
            Kkp(i) = cel(k,  1.0, 1.0, 1.0);
            Ekp(i) = cel(k,  1.0, 1.0, k * k);
        });
        MATAR_FENCE();
        Kk.update_host();
        Ek.update_host();
        Kkp.update_host();
        Ekp.update_host();

        std::printf("Test 2: Legendre relation  E*K' + E'*K - K*K' = pi/2\n");
        std::printf("%6s %12s %12s %12s %12s %14s %12s\n",
                    "k", "K(k)", "E(k)", "K(k')", "E(k')",
                    "Legendre", "Error");

        double max_leg_err = 0.0;
        for (int i = 0; i < NLEG; i++) {
            double leg = Ek.host(i)  * Kkp.host(i)
                       + Ekp.host(i) * Kk.host(i)
                       - Kk.host(i)  * Kkp.host(i);
            double err = std::fabs(leg - PIO2);
            if (err > max_leg_err) max_leg_err = err;
            std::printf("%6.2f %12.7f %12.7f %12.7f %12.7f %14.10f %12.2e\n",
                        test_k[i],
                        Kk.host(i), Ek.host(i),
                        Kkp.host(i), Ekp.host(i),
                        leg, err);
        }

        std::printf("\nMax Legendre relation error: %.2e\n", max_leg_err);

        bool pass = (exact_err < 1e-8) && (max_leg_err < 1e-4);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
