#include <cstdio>
#include <cmath>
#include <matar.h>
#include "eigsrt.hpp"
#include "../jacobi/jacobi.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<double> d(N);
        DFMatrixKokkos<double> v(N, N);

        a.host(1,1) = 2.0;  a.host(1,2) = 1.0;  a.host(1,3) = 0.0;
        a.host(2,1) = 1.0;  a.host(2,2) = 3.0;  a.host(2,3) = 1.0;
        a.host(3,1) = 0.0;  a.host(3,2) = 1.0;  a.host(3,3) = 2.0;

        int nrot = 0;
        jacobi(a, N, d, v, nrot);

        std::printf("Before sorting:\n");
        for (int i = 1; i <= N; i++)
            std::printf("  d(%d) = %16.10f\n", i, d.host(i));

        eigsrt(d, v, N);

        std::printf("\nAfter sorting (descending):\n");
        for (int i = 1; i <= N; i++)
            std::printf("  d(%d) = %16.10f\n", i, d.host(i));

        std::printf("\nSorted eigenvectors (columns of V):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.6f", v.host(i, j));
            std::printf(" ]\n");
        }

        double expected[] = {4.0, 2.0, 1.0};
        double max_err = 0.0;
        for (int i = 0; i < N; i++) {
            double err = std::fabs(d.host(i + 1) - expected[i]);
            if (err > max_err) max_err = err;
        }

        bool order_ok = true;
        for (int i = 1; i < N; i++) {
            if (d.host(i) < d.host(i + 1)) order_ok = false;
        }

        std::printf("\nMax eigenvalue error: %.2e\n", max_err);
        std::printf("Descending order: %s\n", order_ok ? "YES" : "NO");
        std::printf("Test %s\n", (max_err < 1e-10 && order_ok) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
