#include <cstdio>
#include <cmath>
#include <matar.h>
#include "simp1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Simplex LP helper: simp1\n");
        std::printf("========================\n\n");

        // 3x4 tableau (row-major)
        constexpr int mp = 3, np = 4;
        double a[mp * np] = {
             1.0, -2.0,  3.0, -1.0,
             0.5,  4.0, -1.0,  2.0,
            -3.0,  1.0,  5.0,  0.5
        };

        int ll[3] = {0, 1, 2};

        // Find max in row 0 among columns {0,1,2}
        int kp = 0;
        double bmax = 0.0;
        simp1(a, np, 0, ll, 3, 0, kp, bmax);

        std::printf("Row 0: [1.0, -2.0, 3.0, -1.0]\n");
        std::printf("Scan columns {0,1,2}: max at col %d, bmax = %.2f\n", kp, bmax);
        bool pass = (kp == 2) && (std::fabs(bmax - 3.0) < 1.0e-10);

        // Test with absolute value comparison
        simp1(a, np, 0, ll, 3, 1, kp, bmax);
        std::printf("Scan columns {0,1,2} by |val|: max at col %d, |bmax| = %.2f\n",
                    kp, std::fabs(bmax));
        pass = pass && (kp == 2);

        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
