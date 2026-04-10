// General linear least squares (Numerical Recipes LFIT).
//
// Fits y(i) = sum_k a(k)*afunc_k(x(i)) by solving the normal equations
// via Gauss-Jordan elimination. Parameters listed in lista(1..mfit) are
// fitted; the rest are held frozen at their input values.
//
// funcs(x, afunc, ma): evaluates the ma basis functions at x.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "gaussj.hpp"
#include "covsrt.hpp"

using namespace mtr;

template<typename Func>
inline void lfit(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                 DFMatrixKokkos<double>& sig, int ndata,
                 DFMatrixKokkos<double>& a, int ma,
                 DFMatrixKokkos<int>& lista, int mfit,
                 DFMatrixKokkos<double>& covar,
                 double& chisq, Func funcs)
{
    int kk = mfit + 1;
    for (int j = 1; j <= ma; j++) {
        int ihit = 0;
        for (int k = 1; k <= mfit; k++)
            if (lista.host(k) == j) ihit++;
        if (ihit == 0) {
            lista.host(kk) = j;
            kk++;
        } else if (ihit > 1) {
            std::fprintf(stderr, "lfit: improper set in lista\n");
            return;
        }
    }
    if (kk != ma + 1) {
        std::fprintf(stderr, "lfit: improper set in lista\n");
        return;
    }

    double beta_arr[mfit + 1];
    for (int j = 1; j <= mfit; j++) {
        for (int k = 1; k <= mfit; k++)
            covar.host(j, k) = 0.0;
        beta_arr[j] = 0.0;
    }

    DFMatrixKokkos<double> afunc(ma);

    for (int i = 1; i <= ndata; i++) {
        funcs(x.host(i), afunc, ma);
        double ym = y.host(i);
        if (mfit < ma) {
            for (int j = mfit + 1; j <= ma; j++)
                ym -= a.host(lista.host(j)) * afunc.host(lista.host(j));
        }
        double sig2i = 1.0 / (sig.host(i) * sig.host(i));
        for (int j = 1; j <= mfit; j++) {
            double wt = afunc.host(lista.host(j)) * sig2i;
            for (int k = 1; k <= j; k++)
                covar.host(j, k) += wt * afunc.host(lista.host(k));
            beta_arr[j] += ym * wt;
        }
    }

    if (mfit > 1) {
        for (int j = 2; j <= mfit; j++)
            for (int k = 1; k <= j - 1; k++)
                covar.host(k, j) = covar.host(j, k);
    }

    DFMatrixKokkos<double> beta_mat(mfit, 1);
    for (int j = 1; j <= mfit; j++)
        beta_mat.host(j, 1) = beta_arr[j];

    gaussj(covar, mfit, beta_mat, 1);

    for (int j = 1; j <= mfit; j++)
        a.host(lista.host(j)) = beta_mat.host(j, 1);

    chisq = 0.0;
    for (int i = 1; i <= ndata; i++) {
        funcs(x.host(i), afunc, ma);
        double sum = 0.0;
        for (int j = 1; j <= ma; j++)
            sum += a.host(j) * afunc.host(j);
        double r = (y.host(i) - sum) / sig.host(i);
        chisq += r * r;
    }

    covsrt(covar, ma, lista, mfit);
}
