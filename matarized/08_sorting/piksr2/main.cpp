#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <matar.h>
#include "piksr2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 15;

        DFMatrixKokkos<double> arr(N);
        DFMatrixKokkos<double> brr(N);

        std::srand(54321);
        std::printf("Original arrays:\n");
        std::printf("%6s %10s %10s\n", "Index", "arr", "brr");
        for (int i = 1; i <= N; i++) {
            arr.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            brr.host(i) = i * 100.0;
            std::printf("%6d %10.2f %10.2f\n", i, arr.host(i), brr.host(i));
        }

        double orig_arr[N + 1], orig_brr[N + 1];
        for (int i = 1; i <= N; i++) {
            orig_arr[i] = arr.host(i);
            orig_brr[i] = brr.host(i);
        }

        piksr2(N, arr, brr);

        std::printf("\nSorted by arr (insertion sort of 2 arrays):\n");
        std::printf("%6s %10s %10s\n", "Index", "arr", "brr");
        bool sorted = true;
        bool correspondence = true;
        for (int i = 1; i <= N; i++) {
            std::printf("%6d %10.2f %10.2f\n", i, arr.host(i), brr.host(i));
            if (i > 1 && arr.host(i) < arr.host(i - 1)) sorted = false;
        }

        for (int i = 1; i <= N; i++) {
            bool found = false;
            for (int j = 1; j <= N; j++) {
                if (std::fabs(arr.host(i) - orig_arr[j]) < 1e-12 &&
                    std::fabs(brr.host(i) - orig_brr[j]) < 1e-12) {
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
