#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gauleg.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const double pi = 3.14159265358979323846;

        std::printf("Gauss-Legendre quadrature abscissas and weights\n\n");

        // Compute for several orders and test integration accuracy
        int orders[] = {5, 10, 20};
        const int n_orders = 3;

        struct TestFunc {
            const char* name;
            double a, b, exact;
        };

        TestFunc tests[] = {
            {"x^2 on [0,1]",      0.0, 1.0, 1.0 / 3.0},
            {"sin(x) on [0,pi]",  0.0, pi,  2.0},
            {"exp(x) on [0,1]",   0.0, 1.0, std::exp(1.0) - 1.0}
        };

        bool all_pass = true;

        for (int oi = 0; oi < n_orders; oi++) {
            int n = orders[oi];
            DCArrayKokkos<double> x(n);
            DCArrayKokkos<double> w(n);

            gauleg(0.0, 1.0, x, w, n);  // abscissas on [0,1]

            std::printf("--- N = %d ---\n", n);
            std::printf("%4s %18s %18s\n", "i", "x_i", "w_i");
            for (int i = 0; i < n; i++)
                std::printf("%4d %18.14f %18.14f\n", i, x.host(i), w.host(i));

            // Verify weights sum to interval length
            double wsum = 0.0;
            for (int i = 0; i < n; i++)
                wsum += w.host(i);
            std::printf("  Sum of weights = %.14f (should be 1.0)\n\n", wsum);

            // Integration tests using the computed quadrature
            for (int ti = 0; ti < 3; ti++) {
                DCArrayKokkos<double> xq(n), wq(n);
                gauleg(tests[ti].a, tests[ti].b, xq, wq, n);

                double result = 0.0;
                for (int i = 0; i < n; i++) {
                    double xi = xq.host(i);
                    double wi = wq.host(i);
                    double fval;
                    if (ti == 0)      fval = xi * xi;
                    else if (ti == 1) fval = std::sin(xi);
                    else              fval = std::exp(xi);
                    result += wi * fval;
                }
                double err = std::fabs(result - tests[ti].exact);
                std::printf("  N=%2d  %-25s result=%.12f  exact=%.12f  err=%.2e\n",
                            n, tests[ti].name, result, tests[ti].exact, err);
                if (err > 1e-4) all_pass = false;
            }
            std::printf("\n");
        }

        std::printf("Test %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
