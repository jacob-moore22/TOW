#include <cstdio>
#include <cmath>
#include <matar.h>
#include "f1dim.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("1-D line search evaluator (f1dim)\n");
        std::printf("=================================\n\n");

        // 2D Rosenbrock: f(x,y) = (1-x)^2 + 100(y-x^2)^2
        auto rosenbrock = [](const double* v) {
            double a = 1.0 - v[0];
            double b = v[1] - v[0] * v[0];
            return a * a + 100.0 * b * b;
        };

        double p[2]  = {0.0, 0.0};
        double xi[2] = {1.0, 1.0};

        auto f1 = make_f1dim(p, xi, 2, rosenbrock);

        std::printf("Base point p  = (%.1f, %.1f)\n", p[0], p[1]);
        std::printf("Direction  xi = (%.1f, %.1f)\n\n", xi[0], xi[1]);
        std::printf("%10s %18s %18s\n", "t", "f1dim(t)", "rosenbrock(p+t*xi)");

        bool pass = true;
        for (int i = 0; i <= 10; i++) {
            double t = i * 0.2;
            double val = f1(t);
            double pt[2] = { p[0] + t * xi[0], p[1] + t * xi[1] };
            double expected = rosenbrock(pt);
            std::printf("%10.2f %18.10f %18.10f\n", t, val, expected);
            if (std::fabs(val - expected) > 1.0e-12) pass = false;
        }

        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
