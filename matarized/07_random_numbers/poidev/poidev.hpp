#pragma once
#include <cmath>
#include <matar.h>
#include "ran1.hpp"
#include "gammln.hpp"

// Cached state for poidev to avoid recomputation when xm is unchanged.
struct PoidevState {
    double oldm;
    double sq, alxm, g;
    PoidevState() : oldm(-1.0), sq(0.0), alxm(0.0), g(0.0) {}
};

// Poisson-distributed deviate with mean xm, using ran1 as
// the uniform source and gammln for log-gamma evaluation.
inline double poidev(double xm, int& idum, Ran1State& rstate,
                     PoidevState& pstate)
{
    constexpr double PI = 3.141592654;
    double em, t;

    if (xm < 12.0) {
        if (xm != pstate.oldm) {
            pstate.oldm = xm;
            pstate.g = exp(-xm);
        }
        em = -1.0;
        t  =  1.0;
        do {
            em += 1.0;
            t  *= ran1(idum, rstate);
        } while (t > pstate.g);
    } else {
        if (xm != pstate.oldm) {
            pstate.oldm = xm;
            pstate.sq   = sqrt(2.0 * xm);
            pstate.alxm = log(xm);
            pstate.g    = xm * pstate.alxm - gammln(xm + 1.0);
        }
        double y;
        for (;;) {
            y  = tan(PI * ran1(idum, rstate));
            em = pstate.sq * y + xm;
            if (em < 0.0) continue;
            em = static_cast<double>(static_cast<int>(em));
            t  = 0.9 * (1.0 + y * y) *
                 exp(em * pstate.alxm - gammln(em + 1.0) - pstate.g);
            if (ran1(idum, rstate) <= t) break;
        }
    }
    return em;
}
