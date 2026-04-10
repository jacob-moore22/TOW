#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Euler series acceleration for alternating series.
// Call sequentially for jterm = 1, 2, 3, ... passing each term.
// nterm must persist between calls (pass by reference).
inline void eulsum(double& sum, double term, int jterm,
                   DFMatrixKokkos<double>& wksp, int& nterm)
{
    if (jterm == 1) {
        nterm    = 1;
        wksp(1)  = term;
        sum      = 0.5 * term;
    } else {
        double tmp = wksp(1);
        wksp(1)    = term;
        for (int j = 1; j <= nterm; j++) {
            double dum  = wksp(j + 1);
            wksp(j + 1) = 0.5 * (wksp(j) + tmp);
            tmp         = dum;
        }
        if (std::fabs(wksp(nterm + 1)) <= std::fabs(wksp(nterm))) {
            sum += 0.5 * wksp(nterm + 1);
            nterm++;
        } else {
            sum += wksp(nterm + 1);
        }
    }
}
