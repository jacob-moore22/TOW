#include <cstdio>
#include <cmath>
#include <matar.h>
#include "tridag.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Tridiagonal system with known solution u = [1, 2, 3, 4, 5]
        // Diagonals: a(sub), b(main), c(super)
        //  [ 2 -1  0  0  0 ] [1]   [ 0]
        //  [-1  2 -1  0  0 ] [2]   [ 0]
        //  [ 0 -1  2 -1  0 ] [3] = [ 0]
        //  [ 0  0 -1  2 -1 ] [4]   [ 0]
        //  [ 0  0  0 -1  2 ] [5]   [ 6]
        constexpr int N = 5;

        DFMatrixKokkos<double> a(N), b(N), c(N), r(N), u(N);

        double x_exact[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

        for (int i = 1; i <= N; i++) {
            a.host(i) = -1.0;
            b.host(i) =  2.0;
            c.host(i) = -1.0;
        }

        // Compute RHS: r = T * x_exact
        for (int i = 1; i <= N; i++) {
            r.host(i) = b.host(i) * x_exact[i-1];
            if (i > 1) r.host(i) += a.host(i) * x_exact[i-2];
            if (i < N) r.host(i) += c.host(i) * x_exact[i];
        }

        tridag(a, b, c, r, u, N);

        std::printf("Tridiagonal Solver (TRIDAG)\n");
        std::printf("===========================\n\n");

        std::printf("Solution:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(u.host(i) - x_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  u(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, u.host(i), x_exact[i-1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
