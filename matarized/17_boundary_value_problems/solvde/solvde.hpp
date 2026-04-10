#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Forward declarations of helper routines
inline void bksub(int ne, int nb, int jf, int k1, int k2,
                  DFMatrixKokkos<double>& c, int nci, int ncj, int nck);

inline void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
                  DFMatrixKokkos<double>& c, int nci, int ncj, int nck,
                  DFMatrixKokkos<double>& s, int nsi, int nsj);

inline void red(int iz1, int iz2, int jz1, int jz2,
                int jm1, int jm2, int jmf,
                int ic1, int jc1, int jcf, int kc,
                DFMatrixKokkos<double>& c, int nci, int ncj, int nck,
                DFMatrixKokkos<double>& s, int nsi, int nsj);

// General BVP relaxation solver (Numerical Recipes SOLVDE).
// Difeq is a callable with signature:
//   void(int k, int k1, int k2, int jsf,
//        int is1, int isf, int* indexv, int ne,
//        DFMatrixKokkos<double>& s, int nsi, int nsj,
//        DFMatrixKokkos<double>& y, int nyj, int nyk)
template <typename Difeq>
inline void solvde(int itmax, double conv, double slowc,
                   double* scalv, int* indexv,
                   int ne, int nb, int m,
                   DFMatrixKokkos<double>& y, int nyj, int nyk,
                   DFMatrixKokkos<double>& c, int nci, int ncj, int nck,
                   DFMatrixKokkos<double>& s, int nsi, int nsj,
                   Difeq difeq)
{
    constexpr int NMAX = 10;
    double ermax_arr[NMAX];
    int    kmax_arr[NMAX];

    int k1 = 1;
    int k2 = m;
    int nvars = ne * m;

    int j1 = 1;
    int j2 = nb;
    int j3 = nb + 1;
    int j4 = ne;
    int j5 = j4 + j1;
    int j6 = j4 + j2;
    int j7 = j4 + j3;
    int j8 = j4 + j4;
    int j9 = j8 + j1;

    int ic1 = 1;
    int ic2 = ne - nb;
    int ic3 = ic2 + 1;
    int ic4 = ne;

    int jc1 = 1;
    int jcf = ic3;

    for (int it = 1; it <= itmax; it++) {
        int k = k1;
        difeq(k, k1, k2, j9, ic3, ic4, indexv, ne, s, nsi, nsj, y, nyj, nyk);
        pinvs(ic3, ic4, j5, j9, jc1, k1, c, nci, ncj, nck, s, nsi, nsj);

        for (k = k1 + 1; k <= k2; k++) {
            int kp = k - 1;
            difeq(k, k1, k2, j9, ic1, ic4, indexv, ne, s, nsi, nsj, y, nyj, nyk);
            red(ic1, ic4, j1, j2, j3, j4, j9, ic3, jc1, jcf, kp,
                c, nci, ncj, nck, s, nsi, nsj);
            pinvs(ic1, ic4, j3, j9, jc1, k, c, nci, ncj, nck, s, nsi, nsj);
        }

        k = k2 + 1;
        difeq(k, k1, k2, j9, ic1, ic2, indexv, ne, s, nsi, nsj, y, nyj, nyk);
        red(ic1, ic2, j5, j6, j7, j8, j9, ic3, jc1, jcf, k2,
            c, nci, ncj, nck, s, nsi, nsj);
        pinvs(ic1, ic2, j7, j9, jcf, k2 + 1, c, nci, ncj, nck, s, nsi, nsj);
        bksub(ne, nb, jcf, k1, k2, c, nci, ncj, nck);

        double err = 0.0;
        for (int j = 1; j <= ne; j++) {
            int jv = indexv[j];
            ermax_arr[j] = 0.0;
            double errj = 0.0;
            kmax_arr[j] = 0;
            double vmax = 0.0;
            int km = 0;
            for (k = k1; k <= k2; k++) {
                double vz = std::fabs(c.host(j, 1, k));
                if (vz > vmax) {
                    vmax = vz;
                    km   = k;
                }
                errj += vz;
            }
            err += errj / scalv[jv];
            ermax_arr[j] = (km > 0) ? c.host(j, 1, km) / scalv[jv] : 0.0;
            kmax_arr[j]  = km;
        }
        err /= nvars;

        double fac = slowc / std::fmax(slowc, err);
        for (int jv = 1; jv <= ne; jv++) {
            int j = indexv[jv];
            for (k = k1; k <= k2; k++)
                y.host(j, k) -= fac * c.host(jv, 1, k);
        }

        std::printf("  iter %4d  err=%12.6f  fac=%12.6f\n", it, err, fac);

        if (err < conv) return;
    }
    std::fprintf(stderr, "solvde: ITMAX exceeded\n");
}
