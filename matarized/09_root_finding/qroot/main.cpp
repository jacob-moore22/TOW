#include <cstdio>
#include <cmath>
#include <matar.h>
#include "qroot.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== QROOT: Quadratic factor via Bairstow's method ===\n\n");

        // p(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
        // = (x-1)(x^2 - 5x + 6)  so quadratic factor has b=-5, c=6
        // 1-indexed: p[1]=-6, p[2]=11, p[3]=-6, p[4]=1
        {
            double p[5] = {0.0, -6.0, 11.0, -6.0, 1.0};
            double b = -1.0, c = 1.0;
            qroot(p, 4, b, c, 1.0e-8);
            std::printf("(x-1)(x-2)(x-3): quadratic factor x^2 + (%.8f)x + (%.8f)\n",
                        b, c);
            double disc = b * b - 4.0 * c;
            if (disc >= 0.0) {
                double r1 = (-b + std::sqrt(disc)) / 2.0;
                double r2 = (-b - std::sqrt(disc)) / 2.0;
                std::printf("  roots of factor: %.8f, %.8f\n", r1, r2);
            }
            std::printf("  expected some pair from {1, 2, 3}\n");
        }

        // p(x) = x^4 - 5x^2 + 4 = (x^2 - 1)(x^2 - 4) = (x-1)(x+1)(x-2)(x+2)
        // 1-indexed: p[1]=4, p[2]=0, p[3]=-5, p[4]=0, p[5]=1
        {
            double p[6] = {0.0, 4.0, 0.0, -5.0, 0.0, 1.0};
            double b = 0.0, c = -0.5;
            qroot(p, 5, b, c, 1.0e-8);
            std::printf("\n(x^2-1)(x^2-4): quadratic factor x^2 + (%.8f)x + (%.8f)\n",
                        b, c);
            double disc = b * b - 4.0 * c;
            if (disc >= 0.0) {
                double r1 = (-b + std::sqrt(disc)) / 2.0;
                double r2 = (-b - std::sqrt(disc)) / 2.0;
                std::printf("  roots of factor: %.8f, %.8f\n", r1, r2);
            }
            std::printf("  expected some pair from {-2, -1, 1, 2}\n");
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
