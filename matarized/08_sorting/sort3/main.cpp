#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "sort3.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 12;

        DFMatrixKokkos<double> ra(N), rb(N), rc(N);
        DFMatrixKokkos<double> wksp(N);
        DFMatrixKokkos<int>    iwksp(N);

        std::srand(55555);
        std::printf("Original arrays:\n");
        std::printf("%6s %10s %10s %10s\n", "Index", "ra", "rb", "rc");
        for (int i = 1; i <= N; i++) {
            ra.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            rb.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            rc.host(i) = i * 1.0;
            std::printf("%6d %10.2f %10.2f %10.2f\n",
                        i, ra.host(i), rb.host(i), rc.host(i));
        }

        sort3(N, ra, rb, rc, wksp, iwksp);

        std::printf("\nSorted by ra (index-based sort of 3 arrays):\n");
        std::printf("%6s %10s %10s %10s\n", "Index", "ra", "rb", "rc");
        bool passed = true;
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %10.2f %10.2f %10.2f\n",
                        i, ra.host(i), rb.host(i), rc.host(i));
            if (i > 1 && ra.host(i) < ra.host(i - 1)) passed = false;
        }

        std::printf("\nSort order: %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
