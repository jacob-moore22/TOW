// SOR solver using MATAR host-only types (FMatrix: 1-based, column-major).
// Serial execution -- validates MATAR data structures against the C++ baseline.
// Usage: ./sor [JMAX]   (default: 11)

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <matar.h>

using namespace mtr;

void sor(FMatrix<double>& a, FMatrix<double>& b, FMatrix<double>& c,
         FMatrix<double>& d, FMatrix<double>& e, FMatrix<double>& f,
         FMatrix<double>& u, int jmax, double rjac)
{
    const int    MAXITS = 10000;
    const double EPS    = 1.0e-5;

    double anormf = 0.0;
    for (int j = 2; j <= jmax - 1; j++)
        for (int l = 2; l <= jmax - 1; l++)
            anormf += fabs(f(j, l));

    double omega = 1.0;

    for (int n = 1; n <= MAXITS; n++) {
        double anorm = 0.0;

        for (int j = 2; j <= jmax - 1; j++) {
            for (int l = 2; l <= jmax - 1; l++) {
                if ((j + l) % 2 == n % 2) {
                    double resid = a(j, l) * u(j + 1, l)
                                 + b(j, l) * u(j - 1, l)
                                 + c(j, l) * u(j, l + 1)
                                 + d(j, l) * u(j, l - 1)
                                 + e(j, l) * u(j, l)
                                 - f(j, l);
                    anorm    += fabs(resid);
                    u(j, l)  -= omega * resid / e(j, l);
                }
            }
        }

        if (n == 1)
            omega = 1.0 / (1.0 - 0.5 * rjac * rjac);
        else
            omega = 1.0 / (1.0 - 0.25 * rjac * rjac * omega);

        if (n > 1 && anorm < EPS * anormf) {
            printf(" Converged in %d iterations\n", n);
            return;
        }
    }

    fprintf(stderr, "MAXITS exceeded\n");
}

int main(int argc, char* argv[])
{
    const int    JMAX = (argc > 1) ? atoi(argv[1]) : 11;
    const double PI   = 3.141592653589793;

    if (JMAX < 3) {
        fprintf(stderr, "JMAX must be >= 3\n");
        return 1;
    }

    FMatrix<double> a(JMAX, JMAX), b(JMAX, JMAX), c(JMAX, JMAX),
                    d(JMAX, JMAX), e(JMAX, JMAX), f(JMAX, JMAX),
                    u(JMAX, JMAX);

    for (int i = 1; i <= JMAX; i++) {
        for (int j = 1; j <= JMAX; j++) {
            a(i, j) =  1.0;
            b(i, j) =  1.0;
            c(i, j) =  1.0;
            d(i, j) =  1.0;
            e(i, j) = -4.0;
            f(i, j) =  0.0;
            u(i, j) =  0.0;
        }
    }

    int midl = JMAX / 2 + 1;
    f(midl, midl) = 2.0;

    double rjac = cos(PI / JMAX);

    printf(" Grid: %d x %d\n", JMAX, JMAX);
    auto t0 = std::chrono::high_resolution_clock::now();
    sor(a, b, c, d, e, f, u, JMAX, rjac);
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    printf(" Solve time: %.6f s\n", elapsed);
    printf(" Center value u(%d,%d) = %.10f\n", midl, midl, u(midl, midl));

    if (JMAX <= 21) {
        printf("\n SOR Solution:\n");
        for (int i = 1; i <= JMAX; i++) {
            printf(" ");
            for (int j = 1; j <= JMAX; j++)
                printf("%6.2f", u(i, j));
            printf("\n");
        }
    }

    double max_err = 0.0;
    int max_i = 0, max_j = 0;
    for (int i = 2; i <= JMAX - 1; i++) {
        for (int j = 2; j <= JMAX - 1; j++) {
            double val = u(i + 1, j) + u(i - 1, j)
                       + u(i, j + 1) + u(i, j - 1)
                       - 4.0 * u(i, j)
                       - f(i, j);
            if (fabs(val) > fabs(max_err)) {
                max_err = val;
                max_i = i; max_j = j;
            }
        }
    }
    printf(" Max residual error: %.2e at (%d,%d)\n", max_err, max_i, max_j);

    if (JMAX <= 21) {
        printf("\n Test that solution satisfies Difference Eqns:\n");
        for (int i = 2; i <= JMAX - 1; i++) {
            printf("       ");
            for (int j = 2; j <= JMAX - 1; j++) {
                double val = u(i + 1, j) + u(i - 1, j)
                           + u(i, j + 1) + u(i, j - 1)
                           - 4.0 * u(i, j);
                printf("%6.2f", val);
            }
            printf("\n");
        }
    }

    return 0;
}
