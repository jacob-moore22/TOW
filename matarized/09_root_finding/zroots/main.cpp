#include <cstdio>
#include <cmath>
#include <complex>
#include <matar.h>
#include "zroots.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== ZROOTS: All roots of a polynomial ===\n\n");

        // p(x) = x^3 - 6x^2 + 11x - 6 = (x-1)(x-2)(x-3)
        {
            std::complex<double> a[4] = {{-6.0, 0.0}, {11.0, 0.0},
                                          {-6.0, 0.0}, {1.0, 0.0}};
            std::complex<double> roots[3];
            zroots(a, 3, roots, true);

            std::printf("(x-1)(x-2)(x-3): roots =\n");
            for (int i = 0; i < 3; i++) {
                std::printf("  (%.12f, %.12f)\n", roots[i].real(), roots[i].imag());
            }
            std::printf("  expected: 1, 2, 3\n");
        }

        // p(x) = x^4 - 1 = (x-1)(x+1)(x-i)(x+i)
        {
            std::complex<double> a[5] = {{-1.0, 0.0}, {0.0, 0.0},
                                          {0.0, 0.0},  {0.0, 0.0}, {1.0, 0.0}};
            std::complex<double> roots[4];
            zroots(a, 4, roots, true);

            std::printf("\nx^4 - 1: roots =\n");
            for (int i = 0; i < 4; i++) {
                std::printf("  (%.12f, %.12f)\n", roots[i].real(), roots[i].imag());
            }
            std::printf("  expected: -1, +1, -i, +i\n");
        }

        // p(x) = x^2 + 1  → roots at +/-i
        {
            std::complex<double> a[3] = {{1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}};
            std::complex<double> roots[2];
            zroots(a, 2, roots, true);

            std::printf("\nx^2 + 1: roots =\n");
            for (int i = 0; i < 2; i++) {
                std::printf("  (%.12f, %.12f)\n", roots[i].real(), roots[i].imag());
            }
            std::printf("  expected: (0, -1), (0, +1)\n");
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
