#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "qcksrt.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 30;

        DFMatrixKokkos<double> arr(N);

        std::srand(66666);
        std::printf("Original array:\n");
        for (int i = 1; i <= N; i++) {
            arr.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            std::printf("  arr(%2d) = %8.2f\n", i, arr.host(i));
        }

        qcksrt(N, arr);

        std::printf("\nSorted array (quicksort):\n");
        bool passed = true;
        for (int i = 1; i <= N; i++) {
            std::printf("  arr(%2d) = %8.2f\n", i, arr.host(i));
            if (i > 1 && arr.host(i) < arr.host(i - 1)) passed = false;
        }

        std::printf("\nTest %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
