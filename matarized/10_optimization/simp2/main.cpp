#include <cstdio>
#include <cmath>
#include <matar.h>
#include "simp2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Simplex LP helper: simp2 (pivot row selection)\n");
        std::printf("===============================================\n\n");

        // Small tableau: 4 rows x 3 cols (row 0 is objective)
        constexpr int np = 3;
        double a[] = {
            0.0,   0.0,   0.0,   // row 0: objective (unused here)
            6.0,  -2.0,  -1.0,   // row 1: ratio = 6/2 = 3
            4.0,  -1.0,  -3.0,   // row 2: ratio = 4/1 = 4
            8.0,  -4.0,  -2.0    // row 3: ratio = 8/4 = 2  ← pivot
        };

        int l2[] = {1, 2, 3};
        int ip = 0;
        double q1 = 0.0;
        simp2(a, 3, 2, np, l2, 3, ip, 1, q1);

        std::printf("Pivot column kp = 1\n");
        std::printf("Ratios: row1=3.0, row2=4.0, row3=2.0\n");
        std::printf("Selected pivot row ip = %d  (expected 3)\n", ip);
        std::printf("Minimum ratio q1 = %.2f  (expected 2.0)\n\n", q1);

        bool pass = (ip == 3) && (std::fabs(q1 - 2.0) < 1.0e-10);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
