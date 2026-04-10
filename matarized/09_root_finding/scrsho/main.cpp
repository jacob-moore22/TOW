#include <cstdio>
#include <cmath>
#include <matar.h>
#include "scrsho.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== SCRSHO: Screen plot of function ===\n\n");

        std::printf("--- f(x) = sin(x) on [0, 2*pi] ---\n");
        scrsho([](double x) { return std::sin(x); }, 0.0, 2.0 * M_PI);

        std::printf("\n--- f(x) = x^2 - 2 on [-2, 2] ---\n");
        scrsho([](double x) { return x * x - 2.0; }, -2.0, 2.0);

        std::printf("\n--- f(x) = x^3 - x - 1 on [-1, 2] ---\n");
        scrsho([](double x) { return x * x * x - x - 1.0; }, -1.0, 2.0);

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
