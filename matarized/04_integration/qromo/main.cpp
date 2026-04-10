#include <cstdio>
#include <cmath>
#include <matar.h>
#include "qromo.hpp"
#include "midpnt.hpp"
#include "midinf.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const double pi = 3.14159265358979323846;

        std::printf("Romberg open-formula integration (qromo) tests\n\n");

        // Wrap midpnt/midinf as callable rules for qromo
        auto mp_rule = [](auto f, double a, double b, double& s, int n, int& it) {
            midpnt(f, a, b, s, n, it);
        };
        auto mi_rule = [](auto f, double a, double b, double& s, int n, int& it) {
            midinf(f, a, b, s, n, it);
        };

        // --- Test 1: finite integrals with midpnt rule ---
        std::printf("Using midpnt rule:\n");
        std::printf("%-30s %18s %18s %14s\n",
                    "Integral", "Computed", "Exact", "Error");

        auto f_x2  = [](double x) { return x * x; };
        auto f_sin = [](double x) { return std::sin(x); };
        auto f_exp = [](double x) { return std::exp(x); };

        double r1 = qromo(f_x2,  0.0, 1.0, mp_rule);
        double r2 = qromo(f_sin, 0.0, pi,  mp_rule);
        double r3 = qromo(f_exp, 0.0, 1.0, mp_rule);

        double e1 = 1.0 / 3.0;
        double e2 = 2.0;
        double e3 = std::exp(1.0) - 1.0;

        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "x^2 on [0,1]", r1, e1, std::fabs(r1 - e1));
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "sin(x) on [0,pi]", r2, e2, std::fabs(r2 - e2));
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "exp(x) on [0,1]", r3, e3, std::fabs(r3 - e3));

        // --- Test 2: semi-infinite integral with midinf rule ---
        std::printf("\nUsing midinf rule (semi-infinite):\n");

        auto f_exp_neg = [](double x) { return std::exp(-x); };
        double r4 = qromo(f_exp_neg, 1.0, 1.0e30, mi_rule);
        double e4 = std::exp(-1.0);
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "exp(-x) on [1,inf)", r4, e4, std::fabs(r4 - e4));

        auto f_invx2 = [](double x) { return 1.0 / (x * x); };
        double r5 = qromo(f_invx2, 1.0, 1.0e30, mi_rule);
        double e5 = 1.0;
        std::printf("%-30s %18.12f %18.12f %14.2e\n",
                    "1/x^2 on [1,inf)", r5, e5, std::fabs(r5 - e5));

        bool pass = (std::fabs(r1 - e1) < 1e-6) &&
                    (std::fabs(r2 - e2) < 1e-6) &&
                    (std::fabs(r3 - e3) < 1e-6) &&
                    (std::fabs(r4 - e4) < 1e-6) &&
                    (std::fabs(r5 - e5) < 1e-6);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
