#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "sort.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;

        DFMatrixKokkos<double> ra(N);

        std::srand(11111);
        std::printf("Original array:\n");
        for (int i = 1; i <= N; i++) {
            ra.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            std::printf("  ra(%2d) = %8.2f\n", i, ra.host(i));
        }

        sort(N, ra);

        std::printf("\nSorted array (heapsort):\n");
        bool passed = true;
        for (int i = 1; i <= N; i++) {
            std::printf("  ra(%2d) = %8.2f\n", i, ra.host(i));
            if (i > 1 && ra.host(i) < ra.host(i - 1)) passed = false;
        }

        std::printf("\nTest %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
