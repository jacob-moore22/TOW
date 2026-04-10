#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Polynomial division helper (inline copy for self-contained header).
// Divides u[1..n] by v[1..nv], quotient in q[1..n], remainder in r[1..n].
// Coefficients: index 1 = constant, index n = x^(n-1).
inline void qroot_poldiv(const double* u, int n, const double* v, int nv,
                          double* q, double* r)
{
    for (int j = 1; j <= n; j++) {
        r[j] = u[j];
        q[j] = 0.0;
    }
    for (int k = n - nv; k >= 0; k--) {
        q[k + 1] = r[nv + k] / v[nv];
        for (int j = nv + k - 1; j >= k + 1; j--)
            r[j] -= q[k + 1] * v[j - k];
    }
    r[nv] = 0.0;
}

// Bairstow's method: find a quadratic factor (x^2 + b*x + c) of polynomial
// p[1..n]. On entry b and c are initial guesses; on exit they are the
// refined coefficients. p[1] is constant term, p[n] is leading coefficient.
inline void qroot(const double* p, int n, double& b, double& c, double eps)
{
    constexpr int    NMAX  = 101;
    constexpr int    ITMAX = 20;
    constexpr double TINY  = 1.0e-6;

    double d[4];   // d[1]=c, d[2]=b, d[3]=1
    double q[NMAX], qq[NMAX], rem[NMAX];

    d[3] = 1.0;

    for (int iter = 0; iter < ITMAX; iter++) {
        d[2] = b;
        d[1] = c;

        qroot_poldiv(p, n, d, 3, q, rem);
        double s = rem[1];
        double r = rem[2];

        qroot_poldiv(q, n - 1, d, 3, qq, rem);
        double sc = -rem[1];
        double rc = -rem[2];

        // Shift q by one index: q[i+1] = q[i], q[1] = 0
        for (int i = n - 1; i >= 1; i--)
            q[i + 1] = q[i];
        q[1] = 0.0;

        qroot_poldiv(q, n, d, 3, qq, rem);
        double sb = -rem[1];
        double rb = -rem[2];

        double div  = 1.0 / (sb * rc - sc * rb);
        double delb = (r * sc - s * rc) * div;
        double delc = (-r * sb + s * rb) * div;
        b += delb;
        c += delc;

        if ((std::fabs(delb) <= eps * std::fabs(b) || std::fabs(b) < TINY) &&
            (std::fabs(delc) <= eps * std::fabs(c) || std::fabs(c) < TINY))
            return;
    }

    std::fprintf(stderr, "qroot: too many iterations\n");
}
