// Covariance matrix sort (Numerical Recipes COVSRT).
//
// Rearranges the covariance matrix covar(1..ma, 1..ma) so that
// fitted parameters (listed in lista(1..mfit)) are placed in
// their correct rows/columns, and frozen parameters get zero
// rows/columns.

#pragma once
#include <matar.h>

using namespace mtr;

inline void covsrt(DFMatrixKokkos<double>& covar, int ma,
                   DFMatrixKokkos<int>& lista, int mfit)
{
    for (int j = 1; j <= ma - 1; j++)
        for (int i = j + 1; i <= ma; i++)
            covar.host(i, j) = 0.0;

    for (int i = 1; i <= mfit - 1; i++) {
        for (int j = i + 1; j <= mfit; j++) {
            if (lista.host(j) > lista.host(i))
                covar.host(lista.host(j), lista.host(i)) = covar.host(i, j);
            else
                covar.host(lista.host(i), lista.host(j)) = covar.host(i, j);
        }
    }

    double swap = covar.host(1, 1);
    for (int j = 1; j <= ma; j++) {
        covar.host(1, j) = covar.host(j, j);
        covar.host(j, j) = 0.0;
    }
    covar.host(lista.host(1), lista.host(1)) = swap;
    for (int j = 2; j <= mfit; j++)
        covar.host(lista.host(j), lista.host(j)) = covar.host(1, j);

    for (int j = 2; j <= ma; j++)
        for (int i = 1; i <= j - 1; i++)
            covar.host(i, j) = covar.host(j, i);
}
