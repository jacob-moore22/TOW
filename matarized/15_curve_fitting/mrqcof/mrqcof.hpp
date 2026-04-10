// Build normal equations for Levenberg-Marquardt (Numerical Recipes MRQCOF).
//
// Evaluates the linearized fitting matrix alpha(1..mfit, 1..mfit)
// and gradient vector beta(1..mfit) for the chi-squared merit function.
// Called by mrqmin; normally not called directly.
//
// funcs(x, a, ymod, dyda, na): model function returning value and derivatives.

#pragma once
#include <matar.h>

using namespace mtr;

template<typename Func>
inline void mrqcof(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   DFMatrixKokkos<double>& sig, int ndata,
                   DFMatrixKokkos<double>& a, int ma,
                   DFMatrixKokkos<int>& lista, int mfit,
                   DFMatrixKokkos<double>& alpha,
                   DFMatrixKokkos<double>& beta,
                   double& chisq, Func funcs)
{
    for (int j = 1; j <= mfit; j++) {
        for (int k = 1; k <= j; k++)
            alpha.host(j, k) = 0.0;
        beta.host(j) = 0.0;
    }
    chisq = 0.0;

    DFMatrixKokkos<double> dyda(ma);

    for (int i = 1; i <= ndata; i++) {
        double ymod;
        funcs(x.host(i), a, ymod, dyda, ma);
        double sig2i = 1.0 / (sig.host(i) * sig.host(i));
        double dy = y.host(i) - ymod;

        for (int j = 1; j <= mfit; j++) {
            double wt = dyda.host(lista.host(j)) * sig2i;
            for (int k = 1; k <= j; k++)
                alpha.host(j, k) += wt * dyda.host(lista.host(k));
            beta.host(j) += dy * wt;
        }
        chisq += dy * dy * sig2i;
    }

    for (int j = 2; j <= mfit; j++)
        for (int k = 1; k <= j - 1; k++)
            alpha.host(k, j) = alpha.host(j, k);
}
