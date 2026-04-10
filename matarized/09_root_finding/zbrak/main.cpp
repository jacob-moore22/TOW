#include <cstdio>
#include <cmath>
#include <matar.h>
#include "zbrak.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== ZBRAK: Find all root brackets ===\n\n");

        auto f = [](double x) { return std::sin(x); };

        constexpr int MAX_BRACKETS = 20;
        double xb1[MAX_BRACKETS], xb2[MAX_BRACKETS];
        int nb;

        zbrak(f, 0.0, 10.0, 100, xb1, xb2, MAX_BRACKETS, nb);

        std::printf("f(x) = sin(x) on [0, 10], 100 sub-intervals:\n");
        std::printf("Found %d brackets:\n", nb);
        for (int i = 0; i < nb; i++) {
            std::printf("  bracket %d: [%.6f, %.6f]\n", i + 1, xb1[i], xb2[i]);
        }

        std::printf("\nExpected roots near: pi=%.6f, 2*pi=%.6f, 3*pi=%.6f\n",
                    M_PI, 2.0 * M_PI, 3.0 * M_PI);

        auto g = [](double x) { return x * x - 2.0; };
        zbrak(g, -3.0, 3.0, 60, xb1, xb2, MAX_BRACKETS, nb);
        std::printf("\nf(x) = x^2 - 2 on [-3, 3], 60 sub-intervals:\n");
        std::printf("Found %d brackets:\n", nb);
        for (int i = 0; i < nb; i++) {
            std::printf("  bracket %d: [%.6f, %.6f]\n", i + 1, xb1[i], xb2[i]);
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
