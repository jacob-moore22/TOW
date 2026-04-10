#include <cstdio>
#include <matar.h>
#include "eclazz.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 12;

        DFMatrixKokkos<int> nf(N);

        // Elements are equivalent if they have the same value mod 4.
        // Elements: 1..12, so classes should be {4,8,12}, {1,5,9}, {2,6,10}, {3,7,11}
        DFMatrixKokkos<int> values(N);
        for (int i = 1; i <= N; i++) {
            values.host(i) = i;
        }

        auto equiv = [&](int i, int j) -> bool {
            return (values.host(i) % 4) == (values.host(j) % 4);
        };

        eclazz(nf, N, equiv);

        std::printf("Equivalence classes (eclazz, mod 4):\n");
        std::printf("Elements 1..12, equivalent if same value mod 4\n");
        std::printf("Expected classes: {4,8,12}, {1,5,9}, {2,6,10}, {3,7,11}\n\n");
        std::printf("%6s %6s %8s\n", "Elem", "Class", "Val%%4");
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %6d %8d\n", i, nf.host(i), values.host(i) % 4);
        }

        bool passed = true;
        for (int i = 1; i <= N; i++) {
            for (int j = i + 1; j <= N; j++) {
                bool same_class = (nf.host(i) == nf.host(j));
                bool same_mod   = (values.host(i) % 4) == (values.host(j) % 4);
                if (same_class != same_mod) passed = false;
            }
        }

        std::printf("\nTest %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
