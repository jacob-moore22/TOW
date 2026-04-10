// Levenberg-Marquardt nonlinear least squares (Numerical Recipes MRQMIN).
//
// Iteratively minimizes chi-squared for a nonlinear model.
// Call with alamda < 0 to initialize, then repeatedly with alamda > 0.
// Set alamda = 0 for a final call to get the covariance matrix.
//
// funcs(x, a, y, dyda, na): model function returning value and derivatives.
// lista(1..mfit): indices of parameters to fit (others frozen).
//
// Persistent state is kept in static locals matching the original Fortran.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "gaussj.hpp"
#include "covsrt.hpp"
#include "mrqcof.hpp"

using namespace mtr;

template<typename Func>
inline void mrqmin(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   DFMatrixKokkos<double>& sig, int ndata,
                   DFMatrixKokkos<double>& a, int ma,
                   DFMatrixKokkos<int>& lista, int mfit,
                   DFMatrixKokkos<double>& covar,
                   DFMatrixKokkos<double>& alpha,
                   double& chisq, Func funcs, double& alamda)
{
    constexpr int MMAX = 20;
    static double atry_s[MMAX + 1];
    static double beta_s[MMAX + 1];
    static double da_s[MMAX + 1];
    static double ochisq;

    if (alamda < 0.0) {
        int kk = mfit + 1;
        for (int j = 1; j <= ma; j++) {
            int ihit = 0;
            for (int k = 1; k <= mfit; k++)
                if (lista.host(k) == j) ihit++;
            if (ihit == 0) {
                lista.host(kk) = j;
                kk++;
            } else if (ihit > 1) {
                std::fprintf(stderr, "mrqmin: improper permutation in lista\n");
                return;
            }
        }
        if (kk != ma + 1) {
            std::fprintf(stderr, "mrqmin: improper permutation in lista\n");
            return;
        }

        alamda = 0.001;

        DFMatrixKokkos<double> beta_v(mfit);
        mrqcof(x, y, sig, ndata, a, ma, lista, mfit, alpha, beta_v, chisq, funcs);
        for (int j = 1; j <= mfit; j++)
            beta_s[j] = beta_v.host(j);

        ochisq = chisq;
        for (int j = 1; j <= ma; j++)
            atry_s[j] = a.host(j);
    }

    // Build augmented normal equations
    for (int j = 1; j <= mfit; j++) {
        for (int k = 1; k <= mfit; k++)
            covar.host(j, k) = alpha.host(j, k);
        covar.host(j, j) = alpha.host(j, j) * (1.0 + alamda);
        da_s[j] = beta_s[j];
    }

    // Solve: covar * da = beta  (da overwrites beta in the 2D wrapper)
    DFMatrixKokkos<double> da_mat(mfit, 1);
    for (int j = 1; j <= mfit; j++)
        da_mat.host(j, 1) = da_s[j];
    gaussj(covar, mfit, da_mat, 1);
    for (int j = 1; j <= mfit; j++)
        da_s[j] = da_mat.host(j, 1);

    // Final call: sort covariance and return
    if (alamda == 0.0) {
        covsrt(covar, ma, lista, mfit);
        return;
    }

    // Trial step
    for (int j = 1; j <= mfit; j++)
        atry_s[lista.host(j)] = a.host(lista.host(j)) + da_s[j];

    DFMatrixKokkos<double> atry_v(ma);
    for (int j = 1; j <= ma; j++)
        atry_v.host(j) = atry_s[j];

    DFMatrixKokkos<double> da_v(mfit);
    mrqcof(x, y, sig, ndata, atry_v, ma, lista, mfit, covar, da_v, chisq, funcs);

    if (chisq < ochisq) {
        alamda *= 0.1;
        ochisq = chisq;
        for (int j = 1; j <= mfit; j++) {
            for (int k = 1; k <= mfit; k++)
                alpha.host(j, k) = covar.host(j, k);
            beta_s[j] = da_v.host(j);
            a.host(lista.host(j)) = atry_s[lista.host(j)];
        }
    } else {
        alamda *= 10.0;
        chisq = ochisq;
    }
}
