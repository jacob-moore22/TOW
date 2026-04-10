#pragma once
#include <cmath>
#include <matar.h>
#include "ran1.hpp"

// State for gasdev: caches one of the two Box-Muller variates.
struct GasdevState {
    int iset;
    double gset;
    GasdevState() : iset(0), gset(0.0) {}
};

// Gaussian deviate via Box-Muller transform, using ran1 as the
// uniform source.  Returns N(0,1) variates.
inline double gasdev(int& idum, Ran1State& rstate, GasdevState& gstate)
{
    if (gstate.iset == 0) {
        double v1, v2, rsq;
        do {
            v1  = 2.0 * ran1(idum, rstate) - 1.0;
            v2  = 2.0 * ran1(idum, rstate) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        double fac = sqrt(-2.0 * log(rsq) / rsq);
        gstate.gset = v1 * fac;
        gstate.iset = 1;
        return v2 * fac;
    } else {
        gstate.iset = 0;
        return gstate.gset;
    }
}
