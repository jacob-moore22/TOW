#pragma once
#include <cmath>
#include <matar.h>
#include "../mnbrak/mnbrak.hpp"
#include "../brent/brent.hpp"

using namespace mtr;

// Line minimization. Given an n-dimensional point p and direction xi, moves p
// to the minimum of func along the direction xi. On return, xi is rescaled to
// the actual displacement, and fret is the function value at the new point.
//
// Replaces the Fortran COMMON-block / F1DIM pattern by constructing a lambda
// that captures p, xi, n, and func for the 1-D evaluations passed to mnbrak
// and brent.
template <typename Func>
inline void linmin(double* p, double* xi, int n, double& fret, Func func)
{
    constexpr int    NMAX = 50;
    constexpr double TOL  = 1.0e-4;

    // Capture current p and xi for the 1-D line function
    double pcom[NMAX], xicom[NMAX];
    for (int j = 0; j < n; j++) {
        pcom[j]  = p[j];
        xicom[j] = xi[j];
    }

    auto f1dim = [&](double x) -> double {
        double xt[NMAX];
        for (int j = 0; j < n; j++)
            xt[j] = pcom[j] + x * xicom[j];
        return func(xt);
    };

    double ax = 0.0, xx = 1.0, bx = 2.0;
    double fa, fx, fb;
    mnbrak(ax, xx, bx, fa, fx, fb, f1dim);

    double xmin = 0.0;
    fret = brent(ax, xx, bx, f1dim, TOL, xmin);

    for (int j = 0; j < n; j++) {
        xi[j] = xmin * xi[j];
        p[j]  = p[j] + xi[j];
    }
}
