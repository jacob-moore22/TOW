#include <cstdio>
#include <cmath>
#include <matar.h>
#include "splie2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int M = 5, N = 5;
        DFMatrixKokkos<double> x1a(M), x2a(N);
        DFMatrixKokkos<double> ya(M, N), y2a(M, N);

        for (int i = 1; i <= M; i++)
            x1a.host(i) = (i - 1) * M_PI / (M - 1);
        for (int j = 1; j <= N; j++)
            x2a.host(j) = (j - 1) * M_PI / (N - 1);

        for (int i = 1; i <= M; i++)
            for (int j = 1; j <= N; j++)
                ya.host(i, j) = std::sin(x1a.host(i)) * std::cos(x2a.host(j));

        splie2(x1a, x2a, ya, M, N, y2a);

        std::printf("2D spline setup for f(x1,x2) = sin(x1)*cos(x2)\n\n");
        std::printf("Second derivatives y2a (along x2 direction):\n");
        std::printf("%6s ", "i\\j");
        for (int j = 1; j <= N; j++) std::printf("%12d ", j);
        std::printf("\n");
        for (int i = 1; i <= M; i++) {
            std::printf("%6d ", i);
            for (int j = 1; j <= N; j++) {
                std::printf("%12.6f ", y2a.host(i, j));
            }
            std::printf("\n");
        }

        // Verify the spline coefficients are structurally sound:
        // row 1 (sin(0)=0) should have all zero y2a, middle row should be nonzero
        bool ok = true;
        for (int j = 1; j <= N; j++) {
            if (std::fabs(y2a.host(1, j)) > 1e-12) ok = false;
        }

        int imid = (M + 1) / 2;
        bool has_nonzero = false;
        for (int j = 1; j <= N; j++) {
            if (std::fabs(y2a.host(imid, j)) > 1e-6) has_nonzero = true;
        }
        if (!has_nonzero) ok = false;

        std::printf("\nRow 1 all zeros: %s\n", ok ? "yes" : "no");
        std::printf("Row %d has nonzero entries: %s\n", imid, has_nonzero ? "yes" : "no");
        std::printf("Test %s\n", ok ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
