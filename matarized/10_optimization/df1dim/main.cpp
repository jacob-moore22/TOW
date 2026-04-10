#include <cstdio>
#include <cmath>
#include <matar.h>
#include "df1dim.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("1-D directional derivative (df1dim)\n");
        std::printf("====================================\n\n");

        // f(x,y) = x^2 + y^2, grad = (2x, 2y)
        auto grad_sphere = [](const double* v, double* df) {
            df[0] = 2.0 * v[0];
            df[1] = 2.0 * v[1];
        };

        double p[2]  = {1.0, 2.0};
        double xi[2] = {1.0, 0.0};

        auto df1 = make_df1dim(p, xi, 2, grad_sphere);

        std::printf("Base point p  = (%.1f, %.1f)\n", p[0], p[1]);
        std::printf("Direction  xi = (%.1f, %.1f)\n\n", xi[0], xi[1]);
        std::printf("Directional derivative d/dt f(p + t*xi):\n");
        std::printf("  at p: grad f = (2*1, 2*2) = (2, 4)\n");
        std::printf("  dot(grad f, xi) = 2*1 + 4*0 = 2\n\n");

        std::printf("%10s %18s %18s\n", "t", "df1dim(t)", "expected");

        bool pass = true;
        for (int i = 0; i <= 5; i++) {
            double t = i * 0.5;
            double val = df1(t);
            // grad f at (1+t, 2) = (2(1+t), 4), dot with (1,0) = 2(1+t)
            double expected = 2.0 * (1.0 + t);
            std::printf("%10.2f %18.10f %18.10f\n", t, val, expected);
            if (std::fabs(val - expected) > 1.0e-12) pass = false;
        }

        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
