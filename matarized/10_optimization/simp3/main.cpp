#include <cstdio>
#include <cmath>
#include <matar.h>
#include "simp3.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Simplex LP helper: simp3 (pivot operation)\n");
        std::printf("==========================================\n\n");

        // 3x3 tableau
        constexpr int np = 3;
        double a[3 * np] = {
            1.0,  2.0,  3.0,
            4.0,  5.0,  6.0,
            7.0,  8.0, 10.0
        };

        std::printf("Before pivot (ip=1, kp=1):\n");
        for (int i = 0; i < 3; i++) {
            std::printf("  [");
            for (int j = 0; j < 3; j++)
                std::printf(" %8.4f", a[i * np + j]);
            std::printf(" ]\n");
        }

        simp3(a, np, 2, 2, 1, 1);

        std::printf("\nAfter pivot:\n");
        for (int i = 0; i < 3; i++) {
            std::printf("  [");
            for (int j = 0; j < 3; j++)
                std::printf(" %8.4f", a[i * np + j]);
            std::printf(" ]\n");
        }

        // After pivoting on a[1][1]=5:
        // a[1][1] should become 1/5 = 0.2
        bool pass = std::fabs(a[1 * np + 1] - 0.2) < 1.0e-10;
        std::printf("\na[1][1] = %.4f  (expected 0.2)\n", a[1 * np + 1]);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
