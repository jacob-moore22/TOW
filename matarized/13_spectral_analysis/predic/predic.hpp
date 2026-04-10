// Linear prediction using LP coefficients from memcof / fixrts.
//
// data:   input signal of ndata elements (1-based DFMatrixKokkos).
// ndata:  signal length.
// d:      LP coefficients, npoles elements (1-based DFMatrixKokkos).
// npoles: number of LP poles.
// future: output predicted values, nfut elements (1-based DFMatrixKokkos).
// nfut:   number of future values to predict.
//
// Each predicted value depends on the previous, so the outer loop
// is inherently serial.  Runs as a single-thread device kernel (RUN).

#pragma once
#include <matar.h>

using namespace mtr;

inline void predic(DFMatrixKokkos<double>& data, int ndata,
                   DFMatrixKokkos<double>& d, int npoles,
                   DFMatrixKokkos<double>& future, int nfut)
{
    DFMatrixKokkos<double> reg(npoles);

    // Fill shift register from the tail of the data
    DO_ALL(j, 1, npoles, {
        reg(j) = data(ndata + 1 - j);
    });
    MATAR_FENCE();

    // Sequential prediction
    RUN({
        for (int j = 1; j <= nfut; j++) {
            double sum = 0.0;
            for (int k = 1; k <= npoles; k++)
                sum += d(k) * reg(k);
            for (int k = npoles; k >= 2; k--)
                reg(k) = reg(k - 1);
            reg(1) = sum;
            future(j) = sum;
        }
    });
    MATAR_FENCE();
}
