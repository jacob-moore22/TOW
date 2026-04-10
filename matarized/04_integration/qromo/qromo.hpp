#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "polint.hpp"

// Romberg integration on an open interval using a caller-supplied
// quadrature rule.
//
// rule must be callable as rule(func, a, b, s, n, it) — same
// signature as midpnt / midinf.
//
// H(j) decreases by 1/9 (open formula), vs 1/4 for qromb.
template<typename Func, typename Rule>
inline double qromo(Func func, double a, double b, Rule rule,
                    double eps = 1.0e-6, int jmax = 14)
{
    constexpr int K  = 5;
    constexpr int KM = K - 1;
    constexpr int JMAXP = 15;

    double s[JMAXP], h[JMAXP];
    h[0] = 1.0;
    int it = 0;
    double ss = 0.0;

    for (int j = 0; j < jmax; j++) {
        rule(func, a, b, s[j], j + 1, it);
        if (j + 1 >= K) {
            double dss;
            polint(&h[j + 1 - K], &s[j + 1 - K], K, 0.0, ss, dss);
            if (std::fabs(dss) < eps * std::fabs(ss))
                return ss;
        }
        s[j + 1] = s[j];
        h[j + 1] = h[j] / 9.0;
    }
    std::fprintf(stderr, "qromo: too many steps\n");
    return ss;
}
