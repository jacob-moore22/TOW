#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "mdian2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test with odd N
        {
            constexpr int N = 21;
            DFMatrixKokkos<double> x(N);
            double vals[N + 1];

            std::srand(99999);
            std::printf("Test 1: odd N = %d\n", N);
            for (int i = 1; i <= N; i++) {
                x.host(i) = static_cast<double>(std::rand() % 10000) / 100.0;
                vals[i] = x.host(i);
            }

            double xmed;
            mdian2(x, N, xmed);

            std::sort(vals + 1, vals + N + 1);
            double expected = vals[(N + 1) / 2];

            std::printf("Median (mdian2): %.4f\n", xmed);
            std::printf("Expected:        %.4f\n", expected);
            std::printf("Test %s\n\n",
                        std::fabs(xmed - expected) < 1e-6 ? "PASSED" : "FAILED");
        }

        // Test with even N
        {
            constexpr int N = 20;
            DFMatrixKokkos<double> x(N);
            double vals[N + 1];

            std::srand(10101);
            std::printf("Test 2: even N = %d\n", N);
            for (int i = 1; i <= N; i++) {
                x.host(i) = static_cast<double>(std::rand() % 10000) / 100.0;
                vals[i] = x.host(i);
            }

            double xmed;
            mdian2(x, N, xmed);

            std::sort(vals + 1, vals + N + 1);
            double expected = 0.5 * (vals[N / 2] + vals[N / 2 + 1]);

            std::printf("Median (mdian2): %.4f\n", xmed);
            std::printf("Expected:        %.4f\n", expected);
            std::printf("Test %s\n",
                        std::fabs(xmed - expected) < 1e-6 ? "PASSED" : "FAILED");
        }
    }
    MATAR_FINALIZE();
    return 0;
}
