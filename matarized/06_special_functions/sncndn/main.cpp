#include <cstdio>
#include <cmath>
#include <matar.h>
#include "sncndn.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test data from Numerical Recipes FNCVAL.DAT (Jacobian Elliptic Function)
        // Columns: parameter m, argument u, expected sn(u | m)
        constexpr int NTEST = 20;
        double test_m[NTEST] = {
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.5, 0.5, 0.5, 0.5, 0.5,
            1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0
        };
        double test_u[NTEST] = {
             0.1,  0.2,  0.5,  1.0,  2.0,
             0.1,  0.2,  0.5,  1.0,  2.0,
             0.1,  0.2,  0.5,  1.0,  2.0,
             4.0, -0.2, -0.5, -1.0, -2.0
        };
        double expected_sn[NTEST] = {
            0.099833, 0.19867, 0.47943, 0.84147, 0.90930,
            0.099751, 0.19802, 0.47075, 0.80300, 0.99466,
            0.099668, 0.19738, 0.46212, 0.76159, 0.96403,
            0.99933, -0.19738, -0.46212, -0.76159, -0.96403
        };

        DCArrayKokkos<double> m_arr(NTEST), u_arr(NTEST);
        DCArrayKokkos<double> sn_arr(NTEST), cn_arr(NTEST), dn_arr(NTEST);

        for (int i = 0; i < NTEST; i++) {
            m_arr.host(i)  = test_m[i];
            u_arr.host(i)  = test_u[i];
        }
        m_arr.update_device();
        u_arr.update_device();

        FOR_ALL(i, 0, NTEST, {
            double emmc = 1.0 - m_arr(i);
            sncndn(u_arr(i), emmc, sn_arr(i), cn_arr(i), dn_arr(i));
        });
        MATAR_FENCE();
        sn_arr.update_host();
        cn_arr.update_host();
        dn_arr.update_host();

        std::printf("Jacobian Elliptic Functions (MATAR parallel evaluation):\n\n");
        std::printf("%5s %8s %12s %12s %14s %14s\n",
                    "m", "u", "Expected SN", "Computed SN",
                    "SN^2+CN^2", "m*SN^2+DN^2");

        double max_sn_err  = 0.0;
        double max_id1_err = 0.0;
        double max_id2_err = 0.0;

        for (int i = 0; i < NTEST; i++) {
            double sn = sn_arr.host(i);
            double cn = cn_arr.host(i);
            double dn = dn_arr.host(i);
            double m  = test_m[i];

            double sn_err  = std::fabs(sn - expected_sn[i]);
            double id1     = sn * sn + cn * cn;
            double id2     = m * sn * sn + dn * dn;
            double id1_err = std::fabs(id1 - 1.0);
            double id2_err = std::fabs(id2 - 1.0);

            if (sn_err  > max_sn_err)  max_sn_err  = sn_err;
            if (id1_err > max_id1_err) max_id1_err = id1_err;
            if (id2_err > max_id2_err) max_id2_err = id2_err;

            std::printf("%5.1f %8.2f %12.6f %12.6f %14.10f %14.10f\n",
                        m, test_u[i], expected_sn[i], sn, id1, id2);
        }

        std::printf("\nMax |SN error|:         %.2e\n", max_sn_err);
        std::printf("Max |SN^2+CN^2 - 1|:   %.2e\n", max_id1_err);
        std::printf("Max |m*SN^2+DN^2 - 1|: %.2e\n", max_id2_err);

        bool pass = (max_sn_err < 5e-4) && (max_id1_err < 1e-6) && (max_id2_err < 1e-6);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
