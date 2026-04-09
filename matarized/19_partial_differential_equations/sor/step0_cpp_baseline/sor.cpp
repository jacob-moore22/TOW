// Direct C++ translation of Numerical Recipes SOR (sor.f + sor.dem)
// Plain C++ with dynamically-sized arrays -- no external dependencies.
// Usage: ./sor [JMAX]   (default: 11)

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <chrono>

void sor(double* a, double* b, double* c, double* d, double* e,
         double* f, double* u, int jmax, double rjac)
{
    const int    MAXITS = 10000;
    const double EPS    = 1.0e-5;
    auto idx = [jmax](int r, int col) { return r * jmax + col; };

    double anormf = 0.0;
    for (int j = 1; j <= jmax - 2; j++)
        for (int l = 1; l <= jmax - 2; l++)
            anormf += fabs(f[idx(j, l)]);

    double omega = 1.0;

    for (int n = 1; n <= MAXITS; n++) {
        double anorm = 0.0;

        for (int j = 1; j <= jmax - 2; j++) {
            for (int l = 1; l <= jmax - 2; l++) {
                if ((j + l) % 2 == n % 2) {
                    double resid = a[idx(j, l)] * u[idx(j + 1, l)]
                                 + b[idx(j, l)] * u[idx(j - 1, l)]
                                 + c[idx(j, l)] * u[idx(j, l + 1)]
                                 + d[idx(j, l)] * u[idx(j, l - 1)]
                                 + e[idx(j, l)] * u[idx(j, l)]
                                 - f[idx(j, l)];
                    anorm         += fabs(resid);
                    u[idx(j, l)]  -= omega * resid / e[idx(j, l)];
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

    const size_t N = (size_t)JMAX * JMAX;
    std::vector<double> a(N,  1.0), b(N,  1.0), c(N,  1.0),
                        d(N,  1.0), e(N, -4.0), f(N,  0.0), u(N, 0.0);

    auto idx = [JMAX](int r, int col) { return r * JMAX + col; };

    int midl = JMAX / 2;
    f[idx(midl, midl)] = 2.0;

    double rjac = cos(PI / JMAX);

    printf(" Grid: %d x %d\n", JMAX, JMAX);
    auto t0 = std::chrono::high_resolution_clock::now();
    sor(a.data(), b.data(), c.data(), d.data(), e.data(),
        f.data(), u.data(), JMAX, rjac);
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();
    printf(" Solve time: %.6f s\n", elapsed);
    printf(" Center value u(%d,%d) = %.10f\n", midl, midl, u[idx(midl, midl)]);

    if (JMAX <= 21) {
        printf("\n SOR Solution:\n");
        for (int i = 0; i < JMAX; i++) {
            printf(" ");
            for (int j = 0; j < JMAX; j++)
                printf("%6.2f", u[idx(i, j)]);
            printf("\n");
        }
    }

    double max_err = 0.0;
    int max_i = 0, max_j = 0;
    for (int i = 1; i <= JMAX - 2; i++) {
        for (int j = 1; j <= JMAX - 2; j++) {
            double val = u[idx(i + 1, j)] + u[idx(i - 1, j)]
                       + u[idx(i, j + 1)] + u[idx(i, j - 1)]
                       - 4.0 * u[idx(i, j)]
                       - f[idx(i, j)];
            if (fabs(val) > fabs(max_err)) {
                max_err = val;
                max_i = i; max_j = j;
            }
        }
    }
    printf(" Max residual error: %.2e at (%d,%d)\n", max_err, max_i, max_j);

    if (JMAX <= 21) {
        printf("\n Test that solution satisfies Difference Eqns:\n");
        for (int i = 1; i <= JMAX - 2; i++) {
            printf("       ");
            for (int j = 1; j <= JMAX - 2; j++) {
                double val = u[idx(i + 1, j)] + u[idx(i - 1, j)]
                           + u[idx(i, j + 1)] + u[idx(i, j - 1)]
                           - 4.0 * u[idx(i, j)];
                printf("%6.2f", val);
            }
            printf("\n");
        }
    }

    return 0;
}
