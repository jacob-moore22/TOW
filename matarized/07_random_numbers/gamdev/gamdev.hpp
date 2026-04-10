#pragma once
#include <cmath>
#include <matar.h>
#include "ran1.hpp"

// Gamma-distributed deviate of integer order ia (ia >= 1),
// using ran1 as the uniform source.
inline double gamdev(int ia, int& idum, Ran1State& state)
{
    if (ia < 1) return 0.0;

    double x;
    if (ia < 6) {
        x = 1.0;
        for (int j = 0; j < ia; j++)
            x *= ran1(idum, state);
        x = -log(x);
    } else {
        double am = static_cast<double>(ia - 1);
        double s  = sqrt(2.0 * am + 1.0);
        for (;;) {
            double v1 = 2.0 * ran1(idum, state) - 1.0;
            double v2 = 2.0 * ran1(idum, state) - 1.0;
            if (v1 * v1 + v2 * v2 > 1.0) continue;
            double y = v2 / v1;
            x = s * y + am;
            if (x <= 0.0) continue;
            double e = (1.0 + y * y) * exp(am * log(x / am) - s * y);
            if (ran1(idum, state) <= e) break;
        }
    }
    return x;
}
