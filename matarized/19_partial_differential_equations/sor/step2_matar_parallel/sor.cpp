// Parallel SOR solver using MATAR Dual types + Kokkos backend.
// Exploits red-black coloring: within each color, all stencil updates
// are independent, enabling full 2D parallelism via DO_REDUCE_SUM.
// Usage: ./sor [JMAX] [--kokkos-threads=N]   (default JMAX: 11)

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <matar.h>

using namespace mtr;

void sor(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b,
         DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& d,
         DFMatrixKokkos<double>& e, DFMatrixKokkos<double>& f,
         DFMatrixKokkos<double>& u, int jmax, double rjac)
{
    const int    MAXITS = 10000;
    const double EPS    = 1.0e-5;

    double loc_anormf = 0.0, anormf = 0.0;
    DO_REDUCE_SUM(j, 2, jmax - 1,
                  l, 2, jmax - 1,
                  loc_anormf, {
        loc_anormf += fabs(f(j, l));
    }, anormf);

    double omega = 1.0;

    for (int n = 1; n <= MAXITS; n++) {
        int parity = n % 2;

        double loc_anorm = 0.0, anorm = 0.0;
        DO_REDUCE_SUM(j, 2, jmax - 1,
                      l, 2, jmax - 1,
                      loc_anorm, {
            if ((j + l) % 2 == parity) {
                double resid = a(j, l) * u(j + 1, l)
                             + b(j, l) * u(j - 1, l)
                             + c(j, l) * u(j, l + 1)
                             + d(j, l) * u(j, l - 1)
                             + e(j, l) * u(j, l)
                             - f(j, l);
                loc_anorm += fabs(resid);
                u(j, l)    = u(j, l) - omega * resid / e(j, l);
            }
        }, anorm);

        if (n == 1)
            omega = 1.0 / (1.0 - 0.5 * rjac * rjac);
        else
            omega = 1.0 / (1.0 - 0.25 * rjac * rjac * omega);

        if (n > 1 && anorm < EPS * anormf) {
            printf(" Converged in %d iterations\n", n);
            return;
        }
    }

    printf("MAXITS exceeded\n");
}

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int JMAX = 11;
        for (int i = 1; i < argc; i++) {
            if (argv[i][0] != '-') {
                JMAX = atoi(argv[i]);
                break;
            }
        }
        const double PI = 3.141592653589793;

        if (JMAX < 3) {
            fprintf(stderr, "JMAX must be >= 3\n");
            MATAR_FINALIZE();
            return 1;
        }

        DFMatrixKokkos<double> a(JMAX, JMAX), b(JMAX, JMAX), c(JMAX, JMAX),
                               d(JMAX, JMAX), e(JMAX, JMAX), f(JMAX, JMAX),
                               u(JMAX, JMAX);

        for (int i = 1; i <= JMAX; i++) {
            for (int j = 1; j <= JMAX; j++) {
                a.host(i, j) =  1.0;
                b.host(i, j) =  1.0;
                c.host(i, j) =  1.0;
                d.host(i, j) =  1.0;
                e.host(i, j) = -4.0;
                f.host(i, j) =  0.0;
                u.host(i, j) =  0.0;
            }
        }

        int midl = JMAX / 2 + 1;
        f.host(midl, midl) = 2.0;

        a.update_device();
        b.update_device();
        c.update_device();
        d.update_device();
        e.update_device();
        f.update_device();
        u.update_device();

        double rjac = cos(PI / JMAX);

        printf(" Grid: %d x %d\n", JMAX, JMAX);
        auto t0 = std::chrono::high_resolution_clock::now();
        sor(a, b, c, d, e, f, u, JMAX, rjac);
        auto t1 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t1 - t0).count();
        printf(" Solve time: %.6f s\n", elapsed);

        MATAR_FENCE();
        u.update_host();

        printf(" Center value u(%d,%d) = %.10f\n", midl, midl, u.host(midl, midl));

        if (JMAX <= 21) {
            printf("\n SOR Solution:\n");
            for (int i = 1; i <= JMAX; i++) {
                printf(" ");
                for (int j = 1; j <= JMAX; j++)
                    printf("%6.2f", u.host(i, j));
                printf("\n");
            }
        }

        f.update_host();
        double max_err = 0.0;
        int max_i = 0, max_j = 0;
        for (int i = 2; i <= JMAX - 1; i++) {
            for (int j = 2; j <= JMAX - 1; j++) {
                double val = u.host(i + 1, j) + u.host(i - 1, j)
                           + u.host(i, j + 1) + u.host(i, j - 1)
                           - 4.0 * u.host(i, j)
                           - f.host(i, j);
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
                    double val = u.host(i + 1, j) + u.host(i - 1, j)
                               + u.host(i, j + 1) + u.host(i, j - 1)
                               - 4.0 * u.host(i, j);
                    printf("%6.2f", val);
                }
                printf("\n");
            }
        }
    }
    MATAR_FINALIZE();
    return 0;
}
