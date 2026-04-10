#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

namespace link_detail {

inline double alen(double x1, double x2, double y1, double y2)
{
    return std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

} // namespace link_detail

inline void revcst(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   DFMatrixKokkos<int>& iorder, int ncity,
                   DFMatrixKokkos<int>& n, double& de)
{
    using link_detail::alen;
    n.host(3) = 1 + (n.host(1) + ncity - 2) % ncity;
    n.host(4) = 1 + n.host(2) % ncity;

    double xx[5], yy[5];
    for (int j = 1; j <= 4; j++) {
        int ii = iorder.host(n.host(j));
        xx[j] = x.host(ii);
        yy[j] = y.host(ii);
    }
    de = -alen(xx[1], xx[3], yy[1], yy[3])
         -alen(xx[2], xx[4], yy[2], yy[4])
         +alen(xx[1], xx[4], yy[1], yy[4])
         +alen(xx[2], xx[3], yy[2], yy[3]);
}

inline void revers(DFMatrixKokkos<int>& iorder, int ncity,
                   DFMatrixKokkos<int>& n)
{
    int nn = (1 + (n.host(2) - n.host(1) + ncity) % ncity) / 2;
    for (int j = 1; j <= nn; j++) {
        int k = 1 + (n.host(1) + j - 2) % ncity;
        int l = 1 + (n.host(2) - j + ncity) % ncity;
        int itmp = iorder.host(k);
        iorder.host(k) = iorder.host(l);
        iorder.host(l) = itmp;
    }
}

inline void trncst(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   DFMatrixKokkos<int>& iorder, int ncity,
                   DFMatrixKokkos<int>& n, double& de)
{
    using link_detail::alen;
    n.host(4) = 1 + n.host(3) % ncity;
    n.host(5) = 1 + (n.host(1) + ncity - 2) % ncity;
    n.host(6) = 1 + n.host(2) % ncity;

    double xx[7], yy[7];
    for (int j = 1; j <= 6; j++) {
        int ii = iorder.host(n.host(j));
        xx[j] = x.host(ii);
        yy[j] = y.host(ii);
    }
    de = -alen(xx[2], xx[6], yy[2], yy[6])
         -alen(xx[1], xx[5], yy[1], yy[5])
         -alen(xx[3], xx[4], yy[3], yy[4])
         +alen(xx[1], xx[3], yy[1], yy[3])
         +alen(xx[2], xx[4], yy[2], yy[4])
         +alen(xx[5], xx[6], yy[5], yy[6]);
}

inline void trnspt(DFMatrixKokkos<int>& iorder, int ncity,
                   DFMatrixKokkos<int>& n)
{
    DFMatrixKokkos<int> jorder(ncity);

    int m1 = 1 + (n.host(2) - n.host(1) + ncity) % ncity;
    int m2 = 1 + (n.host(5) - n.host(4) + ncity) % ncity;
    int m3 = 1 + (n.host(3) - n.host(6) + ncity) % ncity;

    int nn = 1;
    for (int j = 1; j <= m1; j++) {
        int jj = 1 + (j + n.host(1) - 2) % ncity;
        jorder.host(nn) = iorder.host(jj);
        nn++;
    }
    if (m2 > 0) {
        for (int j = 1; j <= m2; j++) {
            int jj = 1 + (j + n.host(4) - 2) % ncity;
            jorder.host(nn) = iorder.host(jj);
            nn++;
        }
    }
    if (m3 > 0) {
        for (int j = 1; j <= m3; j++) {
            int jj = 1 + (j + n.host(6) - 2) % ncity;
            jorder.host(nn) = iorder.host(jj);
            nn++;
        }
    }
    for (int j = 1; j <= ncity; j++)
        iorder.host(j) = jorder.host(j);
}

inline bool metrop(double de, double t, double ran_val)
{
    return (de < 0.0) || (ran_val < std::exp(-de / t));
}
