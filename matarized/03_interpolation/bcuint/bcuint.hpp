// Bicubic interpolation evaluation (Numerical Recipes BCUINT).
// Given function, gradient, and cross-derivative values at the four
// corners of a grid cell [x1l,x1u] x [x2l,x2u], returns interpolated
// value ansy and partial derivatives ansy1, ansy2 at point (x1, x2).

#pragma once
#include <cstdio>
#include <matar.h>
#include "bcucof.hpp"

using namespace mtr;

inline void bcuint(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& y1,
                   DFMatrixKokkos<double>& y2, DFMatrixKokkos<double>& y12,
                   double x1l, double x1u, double x2l, double x2u,
                   double x1, double x2,
                   double& ansy, double& ansy1, double& ansy2)
{
    DFMatrixKokkos<double> c(4, 4);
    bcucof(y, y1, y2, y12, x1u - x1l, x2u - x2l, c);

    if (x1u == x1l || x2u == x2l) {
        std::fprintf(stderr, "bcuint: bad input (zero cell dimension)\n");
        return;
    }
    double t = (x1 - x1l) / (x1u - x1l);
    double u = (x2 - x2l) / (x2u - x2l);

    ansy  = 0.0;
    ansy2 = 0.0;
    ansy1 = 0.0;
    for (int i = 4; i >= 1; i--) {
        ansy  = t * ansy + ((c.host(i, 4) * u + c.host(i, 3)) * u +
                c.host(i, 2)) * u + c.host(i, 1);
        ansy2 = t * ansy2 + (3.0 * c.host(i, 4) * u +
                2.0 * c.host(i, 3)) * u + c.host(i, 2);
        ansy1 = u * ansy1 + (3.0 * c.host(4, i) * t +
                2.0 * c.host(3, i)) * t + c.host(2, i);
    }
    ansy1 = ansy1 / (x1u - x1l);
    ansy2 = ansy2 / (x2u - x2l);
}
