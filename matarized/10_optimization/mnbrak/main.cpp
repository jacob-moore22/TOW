#include <cstdio>
#include <cmath>
#include <matar.h>
#include "mnbrak.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Bracket a minimum (mnbrak)\n");
        std::printf("=========================\n\n");

        // Test 1: f(x) = (x-3)^2, minimum at x=3
        auto f1 = [](double x) { return (x - 3.0) * (x - 3.0); };
        double ax = 0.0, bx = 1.0, cx = 0.0;
        double fa = 0.0, fb = 0.0, fc = 0.0;
        mnbrak(ax, bx, cx, fa, fb, fc, f1);

        std::printf("Test 1: f(x) = (x-3)^2\n");
        std::printf("  Bracket: ax=%.6f  bx=%.6f  cx=%.6f\n", ax, bx, cx);
        std::printf("  Values:  fa=%.6f  fb=%.6f  fc=%.6f\n", fa, fb, fc);
        bool brackets = (fb <= fa) && (fb <= fc);
        std::printf("  fb <= fa && fb <= fc: %s\n\n", brackets ? "yes" : "no");

        // Test 2: f(x) = cos(x), minimum near pi
        auto f2 = [](double x) { return std::cos(x); };
        ax = 2.0; bx = 3.0;
        mnbrak(ax, bx, cx, fa, fb, fc, f2);

        std::printf("Test 2: f(x) = cos(x)\n");
        std::printf("  Bracket: ax=%.6f  bx=%.6f  cx=%.6f\n", ax, bx, cx);
        std::printf("  Values:  fa=%.6f  fb=%.6f  fc=%.6f\n", fa, fb, fc);
        brackets = (fb <= fa) && (fb <= fc);
        std::printf("  fb <= fa && fb <= fc: %s\n\n", brackets ? "yes" : "no");

        std::printf("Test %s\n", brackets ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
