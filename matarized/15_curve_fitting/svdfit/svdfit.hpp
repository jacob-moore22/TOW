// Least squares fit via Singular Value Decomposition (Numerical Recipes SVDFIT).
//
// Fits y(i) = sum_k a(k)*afunc_k(x(i)) using SVD to handle
// ill-conditioned or rank-deficient design matrices.
// Singular values below TOL * max(w) are zeroed.
//
// funcs(x, afunc, ma): evaluates the ma basis functions at x.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "svdcmp.hpp"
#include "svbksb.hpp"

using namespace mtr;

template<typename Func>
inline void svdfit(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   DFMatrixKokkos<double>& sig, int ndata,
                   DFMatrixKokkos<double>& a, int ma,
                   DFMatrixKokkos<double>& u, DFMatrixKokkos<double>& v,
                   DFMatrixKokkos<double>& w,
                   double& chisq, Func funcs)
{
    constexpr double TOL = 1.0e-5;

    DFMatrixKokkos<double> afunc(ma);
    DFMatrixKokkos<double> b(ndata);

    for (int i = 1; i <= ndata; i++) {
        funcs(x.host(i), afunc, ma);
        double tmp = 1.0 / sig.host(i);
        for (int j = 1; j <= ma; j++)
            u.host(i, j) = afunc.host(j) * tmp;
        b.host(i) = y.host(i) * tmp;
    }

    svdcmp(u, ndata, ma, w, v);

    double wmax = 0.0;
    for (int j = 1; j <= ma; j++)
        if (w.host(j) > wmax) wmax = w.host(j);

    double thresh = TOL * wmax;
    for (int j = 1; j <= ma; j++)
        if (w.host(j) < thresh) w.host(j) = 0.0;

    svbksb(u, w, v, ndata, ma, b, a);

    chisq = 0.0;
    for (int i = 1; i <= ndata; i++) {
        funcs(x.host(i), afunc, ma);
        double sum = 0.0;
        for (int j = 1; j <= ma; j++)
            sum += a.host(j) * afunc.host(j);
        double r = (y.host(i) - sum) / sig.host(i);
        chisq += r * r;
    }
}
