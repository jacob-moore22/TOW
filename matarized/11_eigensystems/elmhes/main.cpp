#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "elmhes.hpp"
#include "../hqr/hqr.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 4;

        DFMatrixKokkos<double> a(N, N);
        DFMatrixKokkos<double> a_hess(N, N);
        DFMatrixKokkos<double> wr(N);
        DFMatrixKokkos<double> wi(N);

        // Non-symmetric matrix with known eigenvalues 1, 2, 3, 4
        // Companion matrix for (l-1)(l-2)(l-3)(l-4) = l^4 - 10l^3 + 35l^2 - 50l + 24
        a.host(1,1) = 0.0;  a.host(1,2) = 0.0;  a.host(1,3) = 0.0;  a.host(1,4) = -24.0;
        a.host(2,1) = 1.0;  a.host(2,2) = 0.0;  a.host(2,3) = 0.0;  a.host(2,4) =  50.0;
        a.host(3,1) = 0.0;  a.host(3,2) = 1.0;  a.host(3,3) = 0.0;  a.host(3,4) = -35.0;
        a.host(4,1) = 0.0;  a.host(4,2) = 0.0;  a.host(4,3) = 1.0;  a.host(4,4) =  10.0;

        std::printf("Reduction to Hessenberg form (ELMHES)\n");
        std::printf("Companion matrix for (l-1)(l-2)(l-3)(l-4)\n");
        std::printf("Expected eigenvalues: 1, 2, 3, 4\n\n");

        std::printf("Input matrix:\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a.host(i, j));
            std::printf(" ]\n");
        }

        elmhes(a, N);

        std::printf("\nAfter ELMHES (upper Hessenberg + multipliers below):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a.host(i, j));
            std::printf(" ]\n");
        }

        // Verify Hessenberg structure: the subdiagonal a(i+1,i) and above
        // are the Hessenberg matrix; below are stored multipliers
        bool hess_ok = true;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                a_hess.host(i, j) = a.host(i, j);
            }
        }

        // Zero out the multiplier storage to get pure Hessenberg form
        for (int i = 3; i <= N; i++)
            for (int j = 1; j <= i - 2; j++)
                a_hess.host(i, j) = 0.0;

        std::printf("\nPure Hessenberg form (multipliers zeroed):\n");
        for (int i = 1; i <= N; i++) {
            std::printf("  [");
            for (int j = 1; j <= N; j++)
                std::printf(" %10.4f", a_hess.host(i, j));
            std::printf(" ]\n");
        }

        // Verify eigenvalues via HQR on the Hessenberg form
        hqr(a_hess, N, wr, wi);

        std::printf("\nEigenvalues via HQR:\n");
        for (int i = 1; i <= N; i++) {
            if (std::fabs(wi.host(i)) < 1e-10)
                std::printf("  lambda(%d) = %16.10f\n", i, wr.host(i));
            else
                std::printf("  lambda(%d) = %16.10f + %16.10f i\n",
                            i, wr.host(i), wi.host(i));
        }

        double evals[N];
        for (int i = 0; i < N; i++) evals[i] = wr.host(i + 1);
        std::sort(evals, evals + N);

        double expected[] = {1.0, 2.0, 3.0, 4.0};
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
