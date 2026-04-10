#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bcucof.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test with f(x1,x2) = x1^2 * x2^2 on cell [0,1] x [0,1]
        // Corners: (0,0), (1,0), (1,1), (0,1)
        // f   = {0, 0, 1, 0}
        // f1  = {0, 0, 2, 0}  (df/dx1 = 2*x1*x2^2)
        // f2  = {0, 0, 2, 0}  (df/dx2 = x1^2*2*x2)
        // f12 = {0, 0, 4, 0}  (d2f/dx1dx2 = 4*x1*x2)

        DFMatrixKokkos<double> yv(4), y1v(4), y2v(4), y12v(4);
        DFMatrixKokkos<double> c(4, 4);

        double yvals[]   = {0, 0, 1, 0};
        double y1vals[]  = {0, 0, 2, 0};
        double y2vals[]  = {0, 0, 2, 0};
        double y12vals[] = {0, 0, 4, 0};

        for (int i = 1; i <= 4; i++) {
            yv.host(i)   = yvals[i - 1];
            y1v.host(i)  = y1vals[i - 1];
            y2v.host(i)  = y2vals[i - 1];
            y12v.host(i) = y12vals[i - 1];
        }

        bcucof(yv, y1v, y2v, y12v, 1.0, 1.0, c);

        std::printf("Bicubic coefficients for f(x1,x2) = x1^2 * x2^2 on [0,1]x[0,1]\n\n");
        std::printf("Coefficient matrix C(i,j):\n");
        for (int i = 1; i <= 4; i++) {
            for (int j = 1; j <= 4; j++) {
                std::printf("%8.3f ", c.host(i, j));
            }
            std::printf("\n");
        }

        // For f = x1^2 * x2^2 expressed in bicubic form,
        // the only nonzero coeff should be c(3,3) = 1
        double err = std::fabs(c.host(3, 3) - 1.0);
        bool ok = (err < 1e-10);
        std::printf("\nc(3,3) = %.8f (expected 1.0, err = %.2e)\n", c.host(3, 3), err);
        std::printf("Test %s\n", ok ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
