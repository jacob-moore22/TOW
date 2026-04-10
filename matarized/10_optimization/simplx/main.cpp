#include <cstdio>
#include <cmath>
#include <matar.h>
#include "simplx.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Simplex method for linear programming (simplx)\n");
        std::printf("===============================================\n\n");

        // Example LP:
        //   Maximize z = x1 + x2 + 3*x3 - 0.5*x4
        //   Subject to:
        //     x1 +   2*x3        <= 740   (m1=1)
        //          2*x2     - 7*x4 <= 0     (m1=2)
        //     x2 -   x3 + 2*x4   <= 0.5   (m1=3; note: all <= so m1=3, m2=0, m3=0)
        //     x1 + x2 + x3 + x4   = 9     (m3=1)
        //
        // Actually, let's use a simpler problem:
        //   Maximize z = 5*x1 + 4*x2
        //   Subject to:
        //     6*x1 + 4*x2 <= 24
        //     x1 + 2*x2   <= 6
        //     x1, x2 >= 0
        //
        // Optimal: x1=3, x2=3/2, z=21 ... let's check:
        //   6*3+4*1.5=24 ✓, 3+3=6 ✓, z=15+6=21 ✓
        //
        // Tableau format: (m+2) x (n+1) = 4 x 3
        //   Row 0:     [ 0,  5,  4]       (objective coefficients)
        //   Row 1:     [24, -6, -4]       (note: NR convention negates a_ij)
        //   Row 2:     [ 6, -1, -2]
        //   Row 3:     auxiliary (zeros, used internally)

        constexpr int m = 2, n = 2;
        constexpr int m1 = 2, m2 = 0, m3 = 0;
        constexpr int np = n + 1;
        constexpr int rows = m + 2;

        // NR convention: a[0][0]=0, a[0][1..n]=c_j,
        //                a[i][0]=b_i, a[i][1..n]=-a_ij (negated!)
        double a[rows * np] = {
             0.0,  5.0,  4.0,
            24.0, -6.0, -4.0,
             6.0, -1.0, -2.0,
             0.0,  0.0,  0.0
        };

        int icase = 0;
        int izrov[n], iposv[m];

        simplx(a, m, n, np, m1, m2, m3, icase, izrov, iposv);

        std::printf("icase = %d  (0 = optimal)\n\n", icase);

        if (icase == 0) {
            std::printf("Optimal z = %.6f  (expected 21.0)\n\n", a[0]);
            std::printf("Solution:\n");
            for (int i = 0; i < m; i++) {
                if (iposv[i] >= 1 && iposv[i] <= n)
                    std::printf("  x%d = %.6f\n", iposv[i], a[(i + 1) * np]);
            }
            std::printf("\n");
        }

        bool pass = (icase == 0) && (std::fabs(a[0] - 21.0) < 1.0e-6);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
