#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bcuint.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test with f(x1,x2) = sin(x1)*sin(x2) on cell [0, pi/2] x [0, pi/2]
        // Corners: (0,0), (pi/2,0), (pi/2,pi/2), (0,pi/2)
        double x1l = 0.0,       x1u = M_PI / 2.0;
        double x2l = 0.0,       x2u = M_PI / 2.0;

        double cx1[] = {x1l, x1u, x1u, x1l};
        double cx2[] = {x2l, x2l, x2u, x2u};

        DFMatrixKokkos<double> yv(4), y1v(4), y2v(4), y12v(4);

        for (int i = 0; i < 4; i++) {
            double s1 = std::sin(cx1[i]), c1 = std::cos(cx1[i]);
            double s2 = std::sin(cx2[i]), c2 = std::cos(cx2[i]);
            yv.host(i + 1)   = s1 * s2;
            y1v.host(i + 1)  = c1 * s2;     // df/dx1
            y2v.host(i + 1)  = s1 * c2;     // df/dx2
            y12v.host(i + 1) = c1 * c2;     // d2f/dx1dx2
        }

        double x1_test[] = {0.3, 0.7, 1.0, 1.2};
        double x2_test[] = {0.4, 0.5, 0.8, 1.0};
        constexpr int NTEST = 4;

        std::printf("Bicubic interpolation of f(x1,x2) = sin(x1)*sin(x2)\n");
        std::printf("Cell: [0, pi/2] x [0, pi/2]\n\n");
        std::printf("%8s %8s %14s %14s %14s\n",
                    "x1", "x2", "bcuint", "exact", "error");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double ansy, ansy1, ansy2;
            bcuint(yv, y1v, y2v, y12v,
                   x1l, x1u, x2l, x2u,
                   x1_test[t], x2_test[t],
                   ansy, ansy1, ansy2);
            double exact = std::sin(x1_test[t]) * std::sin(x2_test[t]);
            double err   = std::fabs(ansy - exact);
            if (err > max_err) max_err = err;
            std::printf("%8.4f %8.4f %14.8f %14.8f %14.2e\n",
                        x1_test[t], x2_test[t], ansy, exact, err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 0.02 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
