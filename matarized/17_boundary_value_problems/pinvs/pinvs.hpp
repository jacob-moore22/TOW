#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Partial Gauss-Jordan elimination for solvde (Numerical Recipes PINVS).
// Operates on rows ie1..ie2, columns je1..jsf of s, storing reduced
// columns into c(*,jc1..,k).
inline void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
                  DFMatrixKokkos<double>& c, int /*nci*/, int /*ncj*/, int /*nck*/,
                  DFMatrixKokkos<double>& s, int /*nsi*/, int /*nsj*/)
{
    constexpr int NMAX = 10;
    double pscl[NMAX + 1];
    int    indxr[NMAX + 1];

    int je2 = je1 + ie2 - ie1;
    int js1 = je2 + 1;

    for (int i = ie1; i <= ie2; i++) {
        double big = 0.0;
        for (int j = je1; j <= je2; j++) {
            double tmp = std::fabs(s.host(i, j));
            if (tmp > big) big = tmp;
        }
        if (big == 0.0) {
            std::fprintf(stderr, "pinvs: singular matrix, row all zero\n");
            return;
        }
        pscl[i]  = 1.0 / big;
        indxr[i] = 0;
    }

    for (int id = ie1; id <= ie2; id++) {
        double piv = 0.0;
        int ipiv = ie1, jpiv = je1;
        for (int i = ie1; i <= ie2; i++) {
            if (indxr[i] == 0) {
                double big = 0.0;
                int jp = je1;
                for (int j = je1; j <= je2; j++) {
                    if (std::fabs(s.host(i, j)) > big) {
                        jp  = j;
                        big = std::fabs(s.host(i, j));
                    }
                }
                if (big * pscl[i] > piv) {
                    ipiv = i;
                    jpiv = jp;
                    piv  = big * pscl[i];
                }
            }
        }
        if (s.host(ipiv, jpiv) == 0.0) {
            std::fprintf(stderr, "pinvs: singular matrix\n");
            return;
        }
        indxr[ipiv] = jpiv;
        double pivinv = 1.0 / s.host(ipiv, jpiv);
        for (int j = je1; j <= jsf; j++)
            s.host(ipiv, j) *= pivinv;
        s.host(ipiv, jpiv) = 1.0;

        for (int i = ie1; i <= ie2; i++) {
            if (indxr[i] != jpiv) {
                if (s.host(i, jpiv) != 0.0) {
                    double dum = s.host(i, jpiv);
                    for (int j = je1; j <= jsf; j++)
                        s.host(i, j) -= dum * s.host(ipiv, j);
                    s.host(i, jpiv) = 0.0;
                }
            }
        }
    }

    int jcoff = jc1 - js1;
    int icoff = ie1 - je1;
    for (int i = ie1; i <= ie2; i++) {
        int irow = indxr[i] + icoff;
        for (int j = js1; j <= jsf; j++)
            c.host(irow, j + jcoff, k) = s.host(i, j);
    }
}
