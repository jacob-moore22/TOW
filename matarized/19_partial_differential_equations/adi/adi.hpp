#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Alternating Direction Implicit method (Numerical Recipes ADI).
//
// Solves: a*u(j-1,l) + (b+e)*u(j,l) + c*u(j+1,l) + d*u(j,l-1) + f*u(j,l+1) + g = 0
// on interior points j,l = 2..jmax-1.
//
// k: number of ADI levels (Wachspress parameters)
// alpha, beta: eigenvalue bounds for the splitting
// eps: convergence tolerance
//
// Requires: tridag.hpp
inline void tridag(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b,
                   DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& r,
                   DFMatrixKokkos<double>& u, int n);

inline void adi(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b,
                DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& d,
                DFMatrixKokkos<double>& e, DFMatrixKokkos<double>& f,
                DFMatrixKokkos<double>& g, DFMatrixKokkos<double>& u,
                int jmax, int k, double alpha, double beta, double eps)
{
    const int MAXITS = 100;

    int k1 = k + 1;
    int nr = 1;
    for (int i = 0; i < k; i++) nr *= 2;

    // Wachspress acceleration parameters
    double* alph = new double[k1 + 1];
    double* bet  = new double[k1 + 1];
    alph[1] = alpha;
    bet[1]  = beta;
    for (int j = 1; j <= k; j++) {
        alph[j + 1] = std::sqrt(alph[j] * bet[j]);
        bet[j + 1]  = 0.5 * (alph[j] + bet[j]);
    }

    double** s_arr = new double*[nr + 1];
    for (int i = 0; i <= nr; i++)
        s_arr[i] = new double[k1 + 1];

    s_arr[1][1] = std::sqrt(alph[k1] * bet[k1]);
    for (int j = 1; j <= k; j++) {
        double ab = alph[k1 - j] * bet[k1 - j];
        int lim = 1;
        for (int i = 0; i < j - 1; i++) lim *= 2;
        for (int n = 1; n <= lim; n++) {
            double disc = std::sqrt(s_arr[n][j] * s_arr[n][j] - ab);
            s_arr[2 * n][j + 1]     = s_arr[n][j] + disc;
            s_arr[2 * n - 1][j + 1] = ab / s_arr[2 * n][j + 1];
        }
    }

    double* r_arr = new double[nr + 1];
    for (int n = 1; n <= nr; n++)
        r_arr[n] = s_arr[n][k1];

    // Cleanup s_arr
    for (int i = 0; i <= nr; i++) delete[] s_arr[i];
    delete[] s_arr;

    // Workspace for psi and tridiagonal solves
    DFMatrixKokkos<double> psi(jmax, jmax);
    int nn = jmax - 2;
    DFMatrixKokkos<double> aa(nn), bb(nn), cc(nn), rr(nn), uu(nn);

    // Initialize psi and compute source norm
    double anormg = 0.0;
    for (int j = 2; j <= jmax - 1; j++) {
        for (int l = 2; l <= jmax - 1; l++) {
            anormg += std::fabs(g.host(j, l));
            psi.host(j, l) = -d.host(j, l) * u.host(j, l - 1)
                            + (r_arr[1] - e.host(j, l)) * u.host(j, l)
                            - f.host(j, l) * u.host(j, l + 1);
        }
    }

    int nits = MAXITS / nr;
    for (int kits = 1; kits <= nits; kits++) {
        for (int n = 1; n <= nr; n++) {
            int next = (n == nr) ? 1 : n + 1;
            double rfact = r_arr[n] + r_arr[next];

            // Row sweeps (x-direction tridiagonal solves)
            for (int l = 2; l <= jmax - 1; l++) {
                for (int j = 2; j <= jmax - 1; j++) {
                    aa.host(j - 1) = a.host(j, l);
                    bb.host(j - 1) = b.host(j, l) + r_arr[n];
                    cc.host(j - 1) = c.host(j, l);
                    rr.host(j - 1) = psi.host(j, l) - g.host(j, l);
                }
                tridag(aa, bb, cc, rr, uu, nn);
                for (int j = 2; j <= jmax - 1; j++)
                    psi.host(j, l) = -psi.host(j, l) + 2.0 * r_arr[n] * uu.host(j - 1);
            }

            // Column sweeps (y-direction tridiagonal solves)
            for (int j = 2; j <= jmax - 1; j++) {
                for (int l = 2; l <= jmax - 1; l++) {
                    aa.host(l - 1) = d.host(j, l);
                    bb.host(l - 1) = e.host(j, l) + r_arr[n];
                    cc.host(l - 1) = f.host(j, l);
                    rr.host(l - 1) = psi.host(j, l);
                }
                tridag(aa, bb, cc, rr, uu, nn);
                for (int l = 2; l <= jmax - 1; l++) {
                    u.host(j, l)   = uu.host(l - 1);
                    psi.host(j, l) = -psi.host(j, l) + rfact * uu.host(l - 1);
                }
            }
        }

        // Check convergence
        double anorm = 0.0;
        for (int j = 2; j <= jmax - 1; j++) {
            for (int l = 2; l <= jmax - 1; l++) {
                double resid = a.host(j, l) * u.host(j - 1, l)
                             + (b.host(j, l) + e.host(j, l)) * u.host(j, l)
                             + c.host(j, l) * u.host(j + 1, l)
                             + d.host(j, l) * u.host(j, l - 1)
                             + f.host(j, l) * u.host(j, l + 1)
                             + g.host(j, l);
                anorm += std::fabs(resid);
            }
        }
        if (anorm < eps * anormg) {
            std::printf("  ADI converged in %d cycles (%d sweeps)\n",
                        kits, kits * nr);
            delete[] alph;
            delete[] bet;
            delete[] r_arr;
            return;
        }
    }

    delete[] alph;
    delete[] bet;
    delete[] r_arr;
    std::fprintf(stderr, "adi: MAXITS exceeded\n");
}
