#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "hqr.hpp"
#include "../balanc/balanc.hpp"
#include "../elmhes/elmhes.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 3;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<double> wr(N);
        DFMatrixKokkos<double> wi(N);

        // Companion matrix for (lambda-1)(lambda-2)(lambda-3)
        // = lambda^3 - 6*lambda^2 + 11*lambda - 6
        // Eigenvalues: 1, 2, 3
        a.host(1,1) = 0.0;  a.host(1,2) = 0.0;  a.host(1,3) = 6.0;
        a.host(2,1) = 1.0;  a.host(2,2) = 0.0;  a.host(2,3) = -11.0;
        a.host(3,1) = 0.0;  a.host(3,2) = 1.0;  a.host(3,3) = 6.0;

        std::printf("BALANC + ELMHES + HQR eigenvalue pipeline\n");
        std::printf("Companion matrix for (l-1)(l-2)(l-3)\n");
        std::printf("Expected eigenvalues: 1, 2, 3 (real)\n\n");

        std::printf("Input matrix:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a.host(i, j));
            std::printf(" ]\n");
        }

        balanc(a, N);
        std::printf("\nAfter BALANC:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a.host(i, j));
            std::printf(" ]\n");
        }

        elmhes(a, N);
        std::printf("\nAfter ELMHES (Hessenberg form):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a.host(i, j));
            std::printf(" ]\n");
        }

        hqr(a, N, wr, wi);

        std::printf("\nComputed eigenvalues:\n");
        for (int i = 1; i <= N; i++) {
            if (std::fabs(wi.host(i)) < 1e-10) {
                std::printf("  lambda(%d) = %16.10f\n", i, wr.host(i));
            } else {
                std::printf("  lambda(%d) = %16.10f + %16.10f i\n",
                            i, wr.host(i), wi.host(i));
            }
        }

        // Sort real parts and compare to expected
        double evals[N];
        for (int i = 0; i < N; i++) evals[i] = wr.host(i + 1);
        std::sort(evals, evals + N);

        double expected[] = {1.0, 2.0, 3.0};
        double max_err = 0.0;
        bool all_real = true;
        for (int i = 0; i < N; i++) {
            double err = std::fabs(evals[i] - expected[i]);
            if (err > max_err) max_err = err;
            if (std::fabs(wi.host(i + 1)) > 1e-10) all_real = false;
        }

        std::printf("\nAll eigenvalues real: %s\n", all_real ? "YES" : "NO");
        std::printf("Max eigenvalue error: %.2e\n", max_err);
        std::printf("Test %s\n", (max_err < 1e-8 && all_real) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
