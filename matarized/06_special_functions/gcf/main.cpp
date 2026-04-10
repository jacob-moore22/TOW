#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gcf.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NTEST = 4;
        double test_a[NTEST] = {1.0, 2.0, 3.0, 5.0};
        double test_x[NTEST] = {2.0, 4.0, 5.0, 8.0};

        DCArrayKokkos<double> a_d(NTEST);
        DCArrayKokkos<double> x_d(NTEST);
        DCArrayKokkos<double> result(NTEST);
        DCArrayKokkos<double> gln_out(NTEST);

        for (int i = 0; i < NTEST; i++) {
            a_d.host(i) = test_a[i];
            x_d.host(i) = test_x[i];
        }
        a_d.update_device();
        x_d.update_device();

        FOR_ALL(i, 0, NTEST, {
            double gammcf_val, gln;
            gcf(gammcf_val, a_d(i), x_d(i), gln);
            result(i)  = gammcf_val;
            gln_out(i) = gln;
        });
        MATAR_FENCE();
        result.update_host();
        gln_out.update_host();

        std::printf("Incomplete gamma continued-fraction Q(a,x) via gcf:\n");
        std::printf("%8s %8s %18s %18s\n", "a", "x", "Q(a,x)", "GLN");
        for (int i = 0; i < NTEST; i++) {
            std::printf("%8.2f %8.2f %18.10f %18.10f\n",
                        test_a[i], test_x[i], result.host(i), gln_out.host(i));
        }

        double expected = exp(-2.0);
        double err = fabs(result.host(0) - expected);
        std::printf("\ngcf(1.0, 2.0) = %.10f  expected Q(1,2) = %.10f  err = %.2e\n",
                    result.host(0), expected, err);
        std::printf("Test %s\n", err < 1e-6 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
