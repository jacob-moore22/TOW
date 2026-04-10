#include <cstdio>
#include <cmath>
#include <matar.h>
#include "vander.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Vandermonde system: fit polynomial through 4 points.
        // Points: (1, 10), (2, 26), (3, 58), (4, 112)
        // Polynomial: p(x) = 1 + 2x + 3x^2 + 4x^3  =>  w = [1, 2, 3, 4] (but
        //   Vandermonde ordering maps differently)
        //
        // Actually, the NR VANDER solves: sum x(k)^(i-1) * w(i) = q(k)
        // i.e. w are coefficients of polynomial in the monomial basis:
        //   p(x) = w(1) + w(2)*x + w(3)*x^2 + ...
        constexpr int N = 4;

        DFMatrixKokkos<double> x(N), w(N), q(N);

        double x_vals[] = {1.0, 2.0, 3.0, 4.0};
        // p(x) = 1 + 2x + 3x^2 => p(1)=6, p(2)=17, p(3)=34, p(4)=57
        // Using w_exact = [1, 2, 3, 0] for a degree-2 poly with 4 points
        // Instead use 4 exact coefficients: w = [1, 2, 3, 4]
        // p(1) = 1+2+3+4 = 10
        // p(2) = 1+4+12+32 = 49
        // p(3) = 1+6+27+108 = 142
        // p(4) = 1+8+48+256 = 313
        double w_exact[] = {1.0, 2.0, 3.0, 4.0};
        double q_vals[4];
        for (int k = 0; k < N; k++) {
            double xk = x_vals[k];
            q_vals[k] = 0.0;
            double xpow = 1.0;
            for (int i = 0; i < N; i++) {
                q_vals[k] += w_exact[i] * xpow;
                xpow *= xk;
            }
        }

        for (int i = 1; i <= N; i++) {
            x.host(i) = x_vals[i-1];
            q.host(i) = q_vals[i-1];
        }

        vander(x, w, q, N);

        std::printf("Vandermonde Solver (VANDER)\n");
        std::printf("===========================\n\n");

        std::printf("RHS values q: ");
        for (int i = 0; i < N; i++) std::printf("%.1f ", q_vals[i]);
        std::printf("\n\n");

        std::printf("Polynomial coefficients:\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double err = std::fabs(w.host(i) - w_exact[i-1]);
            if (err > max_err) max_err = err;
            std::printf("  w(%d) = %16.10f  (exact = %.1f, err = %.2e)\n",
                        i, w.host(i), w_exact[i-1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-8 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
