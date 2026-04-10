#include <cstdio>
#include <cmath>
#include <matar.h>
#include "poldiv.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N  = 6;
        constexpr int NV = 4;

        DFMatrixKokkos<double> u(N), v(NV), q(N), r(N);

        // (x-1)^5 = -1 + 5x - 10x^2 + 10x^3 - 5x^4 + x^5
        u(1) = -1.0;  u(2) =  5.0;  u(3) = -10.0;
        u(4) = 10.0;  u(5) = -5.0;  u(6) =   1.0;

        // (x+1)^3 = 1 + 3x + 3x^2 + x^3
        v(1) = 1.0;  v(2) = 3.0;  v(3) = 3.0;  v(4) = 1.0;

        poldiv(u, N, v, NV, q, r);

        std::printf("Polynomial division: (x-1)^5 / (x+1)^3\n\n");

        std::printf("           %10s%10s%10s%10s%10s%10s\n",
                    "x^0", "x^1", "x^2", "x^3", "x^4", "x^5");

        std::printf("Quotient:  ");
        for (int i = 1; i <= N; i++) std::printf("%10.2f", q(i));
        std::printf("\nExpected:  ");
        std::printf("%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",
                    31.0, -8.0, 1.0, 0.0, 0.0, 0.0);

        std::printf("\n           %10s%10s%10s%10s\n",
                    "x^0", "x^1", "x^2", "x^3");
        std::printf("Remainder: ");
        for (int i = 1; i <= NV; i++) std::printf("%10.2f", r(i));
        std::printf("\nExpected:  ");
        std::printf("%10.2f%10.2f%10.2f%10.2f\n",
                    -32.0, -80.0, -80.0, 0.0);

        double q_exp[] = {31.0, -8.0, 1.0, 0.0, 0.0, 0.0};
        double r_exp[] = {-32.0, -80.0, -80.0, 0.0};

        double max_err = 0.0;
        for (int i = 0; i < N; i++) {
            double err = std::fabs(q(i + 1) - q_exp[i]);
            if (err > max_err) max_err = err;
        }
        for (int i = 0; i < NV; i++) {
            double err = std::fabs(r(i + 1) - r_exp[i]);
            if (err > max_err) max_err = err;
        }

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
