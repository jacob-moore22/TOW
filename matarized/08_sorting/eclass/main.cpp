#include <cstdio>
#include <matar.h>
#include "eclass.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 10;
        constexpr int M = 6;

        DFMatrixKokkos<int> nf(N);
        DFMatrixKokkos<int> lista(M);
        DFMatrixKokkos<int> listb(M);

        // {1,2}, {3,4}, {5,6}, {2,3}, {7,8}, {9,10}
        int pairs_a[] = {0, 1, 3, 5, 2, 7, 9};
        int pairs_b[] = {0, 2, 4, 6, 3, 8, 10};
        for (int i = 1; i <= M; i++) {
            lista.host(i) = pairs_a[i];
            listb.host(i) = pairs_b[i];
        }

        eclass(nf, N, lista, listb, M);

        std::printf("Equivalence classes (eclass):\n");
        std::printf("Pairs: (1,2) (3,4) (5,6) (2,3) (7,8) (9,10)\n");
        std::printf("Expected: {1,2,3,4} {5,6} {7,8} {9,10}\n\n");
        std::printf("%6s %6s\n", "Elem", "Class");
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %6d\n", i, nf.host(i));
        }

        // Verify: 1,2,3,4 should share one class; 5,6 another; 7,8 another; 9,10 another
        bool passed = true;
        if (nf.host(1) != nf.host(2) || nf.host(1) != nf.host(3) || nf.host(1) != nf.host(4))
            passed = false;
        if (nf.host(5) != nf.host(6))
            passed = false;
        if (nf.host(7) != nf.host(8))
            passed = false;
        if (nf.host(9) != nf.host(10))
            passed = false;
        if (nf.host(1) == nf.host(5) || nf.host(1) == nf.host(7) || nf.host(1) == nf.host(9))
            passed = false;
        if (nf.host(5) == nf.host(7) || nf.host(5) == nf.host(9))
            passed = false;
        if (nf.host(7) == nf.host(9))
            passed = false;

        std::printf("\nTest %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
