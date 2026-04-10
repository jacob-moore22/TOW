#include <cstdio>
#include <cmath>
#include <matar.h>
#include "probks.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // probks(0) should be 1, probks(large) should be ~0
        double test_vals[] = {0.001, 0.5, 1.0, 1.5, 2.0, 3.0};
        int ntests = 6;

        std::printf("Kolmogorov-Smirnov probability function:\n");
        std::printf("%12s %18s\n", "lambda", "Q_KS(lambda)");

        for (int i = 0; i < ntests; i++) {
            double p = probks(test_vals[i]);
            std::printf("%12.4f %18.10f\n", test_vals[i], p);
        }

        // Known: probks(0.001) ~ 1.0, probks(3.0) ~ 0.0
        bool pass = fabs(probks(0.001) - 1.0) < 0.01 &&
                    probks(3.0) < 1e-6;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
