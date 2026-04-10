#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "mdian1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test with odd N
        {
            constexpr int N = 11;
            DFMatrixKokkos<double> x(N);
            double vals[N + 1];

            std::srand(77777);
            std::printf("Test 1: odd N = %d\n", N);
            for (int i = 1; i <= N; i++) {
                x.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
                vals[i] = x.host(i);
                std::printf("  x(%2d) = %8.2f\n", i, x.host(i));
            }

            double xmed;
            mdian1(x, N, xmed);

            std::sort(vals + 1, vals + N + 1);
            double expected = vals[(N + 1) / 2];

            std::printf("Median (mdian1): %.4f\n", xmed);
            std::printf("Expected:        %.4f\n", expected);
            std::printf("Test %s\n\n", std::fabs(xmed - expected) < 1e-10 ? "PASSED" : "FAILED");
        }

        // Test with even N
        {
            constexpr int N = 10;
            DFMatrixKokkos<double> x(N);
            double vals[N + 1];

            std::srand(88888);
            std::printf("Test 2: even N = %d\n", N);
            for (int i = 1; i <= N; i++) {
                x.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
                vals[i] = x.host(i);
                std::printf("  x(%2d) = %8.2f\n", i, x.host(i));
            }

            double xmed;
            mdian1(x, N, xmed);

            std::sort(vals + 1, vals + N + 1);
            double expected = 0.5 * (vals[N / 2] + vals[N / 2 + 1]);

            std::printf("Median (mdian1): %.4f\n", xmed);
            std::printf("Expected:        %.4f\n", expected);
            std::printf("Test %s\n", std::fabs(xmed - expected) < 1e-10 ? "PASSED" : "FAILED");
        }
    }
    MATAR_FINALIZE();
    return 0;
}
