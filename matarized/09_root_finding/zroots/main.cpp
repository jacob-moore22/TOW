#include <cstdio>
#include <cmath>
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
            const int m = 3;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1) = -6.0; a.host(2) = 0.0;   // coef 0
            a.host(3) = 11.0; a.host(4) = 0.0;   // coef 1
            a.host(5) = -6.0; a.host(6) = 0.0;   // coef 2
            a.host(7) =  1.0; a.host(8) = 0.0;   // coef 3

            DFMatrixKokkos<double> roots(2 * m);
            zroots(a, m, roots, true);

            std::printf("(x-1)(x-2)(x-3): roots =\n");
            for (int i = 0; i < m; i++)
                std::printf("  (%.12f, %.12f)\n",
                            roots.host(2 * i + 1), roots.host(2 * i + 2));
            std::printf("  expected: 1, 2, 3\n");
        }

        // p(x) = x^4 - 1 = (x-1)(x+1)(x-i)(x+i)
        {
            const int m = 4;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1)  = -1.0; a.host(2)  = 0.0;
            a.host(3)  =  0.0; a.host(4)  = 0.0;
            a.host(5)  =  0.0; a.host(6)  = 0.0;
            a.host(7)  =  0.0; a.host(8)  = 0.0;
            a.host(9)  =  1.0; a.host(10) = 0.0;

            DFMatrixKokkos<double> roots(2 * m);
            zroots(a, m, roots, true);

            std::printf("\nx^4 - 1: roots =\n");
            for (int i = 0; i < m; i++)
                std::printf("  (%.12f, %.12f)\n",
                            roots.host(2 * i + 1), roots.host(2 * i + 2));
            std::printf("  expected: -1, +1, -i, +i\n");
        }

        // p(x) = x^2 + 1  → roots at +/-i
        {
            const int m = 2;
            DFMatrixKokkos<double> a(2 * (m + 1));
            a.host(1) = 1.0; a.host(2) = 0.0;
            a.host(3) = 0.0; a.host(4) = 0.0;
            a.host(5) = 1.0; a.host(6) = 0.0;

            DFMatrixKokkos<double> roots(2 * m);
            zroots(a, m, roots, true);

            std::printf("\nx^2 + 1: roots =\n");
            for (int i = 0; i < m; i++)
                std::printf("  (%.12f, %.12f)\n",
                            roots.host(2 * i + 1), roots.host(2 * i + 2));
            std::printf("  expected: (0, -1), (0, +1)\n");
        }

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
