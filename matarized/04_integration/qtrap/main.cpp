#include <cstdio>
#include <cmath>
#include <matar.h>
#include "qtrap.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const double pi = 3.14159265358979323846;

        std::printf("Iterated trapezoidal (qtrap) integration tests\n\n");
        std::printf("%-30s %18s %18s %14s\n",
                    "Integral", "Computed", "Exact", "Error");

        auto f_x2  = [](double x) { return x * x; };
        auto f_sin = [](double x) { return std::sin(x); };
        auto f_exp = [](double x) { return std::exp(x); };

        double r1 = qtrap(f_x2,  0.0, 1.0);
        double r2 = qtrap(f_sin, 0.0, pi);
        double r3 = qtrap(f_exp, 0.0, 1.0);

        double e1 = 1.0 / 3.0;
        double e2 = 2.0;
        double e3 = std::exp(1.0) - 1.0;

        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "x^2 on [0,1]", r1, e1, std::fabs(r1 - e1));
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "sin(x) on [0,pi]", r2, e2, std::fabs(r2 - e2));
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "exp(x) on [0,1]", r3, e3, std::fabs(r3 - e3));

        bool pass = (std::fabs(r1 - e1) < 1e-6) &&
                    (std::fabs(r2 - e2) < 1e-6) &&
                    (std::fabs(r3 - e3) < 1e-6);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
