#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <matar.h>
#include "sort2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 15;

        DFMatrixKokkos<double> ra(N);
        DFMatrixKokkos<double> rb(N);

        std::srand(22222);
        std::printf("Original arrays:\n");
        std::printf("%6s %10s %10s\n", "Index", "ra", "rb");
        for (int i = 1; i <= N; i++) {
            ra.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            rb.host(i) = i * 10.0;
            std::printf("%6d %10.2f %10.2f\n", i, ra.host(i), rb.host(i));
        }

        double orig_ra[N + 1], orig_rb[N + 1];
        for (int i = 1; i <= N; i++) {
            orig_ra[i] = ra.host(i);
            orig_rb[i] = rb.host(i);
        }

        sort2(N, ra, rb);

        std::printf("\nSorted by ra (heapsort of 2 arrays):\n");
        std::printf("%6s %10s %10s\n", "Index", "ra", "rb");
        bool sorted = true;
        bool correspondence = true;
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %10.2f %10.2f\n", i, ra.host(i), rb.host(i));
            if (i > 1 && ra.host(i) < ra.host(i - 1)) sorted = false;
        }

        for (int i = 1; i <= N; i++) {
            bool found = false;
            for (int j = 1; j <= N; j++) {
                if (std::fabs(ra.host(i) - orig_ra[j]) < 1e-12 &&
                    std::fabs(rb.host(i) - orig_rb[j]) < 1e-12) {
                    found = true;
                    break;
                }
            }
            if (!found) correspondence = false;
        }

        std::printf("\nSort order: %s\n", sorted ? "PASSED" : "FAILED");
        std::printf("Pair correspondence: %s\n", correspondence ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
