#pragma once
#include <cmath>
#include <matar.h>
#include "ran1.hpp"
#include "gammln.hpp"

// Cached state for bnldev to avoid recomputation when n or p unchanged.
struct BnldevState {
    int    nold;
    double pold;
    double en, oldg, pc, plog, pclog;
    BnldevState() : nold(-1), pold(-1.0), en(0.0), oldg(0.0),
                     pc(0.0), plog(0.0), pclog(0.0) {}
};

// Binomial deviate: number of successes in n trials with probability pp,
// using ran1 as the uniform source and gammln for log-gamma evaluation.
inline double bnldev(double pp, int n, int& idum, Ran1State& rstate,
                     BnldevState& bstate)
{
    constexpr double PI = 3.141592654;
    double p  = (pp <= 0.5) ? pp : 1.0 - pp;
    double am = n * p;
    double result;

    if (n < 25) {
        result = 0.0;
        for (int j = 0; j < n; j++) {
            if (ran1(idum, rstate) < p) result += 1.0;
        }
    } else if (am < 1.0) {
        double g = exp(-am);
        double t = 1.0;
        int j;
        for (j = 0; j <= n; j++) {
            t *= ran1(idum, rstate);
            if (t < g) break;
        }
        if (j > n) j = n;
        result = static_cast<double>(j);
    } else {
        if (n != bstate.nold) {
            bstate.en   = static_cast<double>(n);
            bstate.oldg = gammln(bstate.en + 1.0);
            bstate.nold = n;
        }
        if (p != bstate.pold) {
            bstate.pc    = 1.0 - p;
            bstate.plog  = log(p);
            bstate.pclog = log(bstate.pc);
            bstate.pold  = p;
        }
        double sq = sqrt(2.0 * am * bstate.pc);
        double em, y, t;
        for (;;) {
            y  = tan(PI * ran1(idum, rstate));
            em = sq * y + am;
            if (em < 0.0 || em >= bstate.en + 1.0) continue;
            em = static_cast<double>(static_cast<int>(em));
            t  = 1.2 * sq * (1.0 + y * y) *
                 exp(bstate.oldg - gammln(em + 1.0) -
                     gammln(bstate.en - em + 1.0) +
                     em * bstate.plog +
                     (bstate.en - em) * bstate.pclog);
            if (ran1(idum, rstate) <= t) break;
        }
        result = em;
    }

    if (p != pp) result = n - result;
    return result;
}
