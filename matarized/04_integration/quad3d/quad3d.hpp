#pragma once
#include <matar.h>
#include "qgaus.hpp"

// Three-dimensional integration using nested 10-point Gauss-Legendre
// quadrature (qgaus).
//
// Computes integral of func(x,y,z) over the region:
//   x  in [x1, x2]
//   y  in [y1(x), y2(x)]
//   z  in [z1(x,y), z2(x,y)]
//
// All function arguments are callables:
//   func : (double, double, double) -> double
//   y1   : (double) -> double
//   y2   : (double) -> double
//   z1   : (double, double) -> double
//   z2   : (double, double) -> double
template<typename Func3D,
         typename Y1Func, typename Y2Func,
         typename Z1Func, typename Z2Func>
inline double quad3d(Func3D func, double x1, double x2,
                     Y1Func y1, Y2Func y2,
                     Z1Func z1, Z2Func z2)
{
    auto h = [&](double x) {
        auto g = [&](double y) {
            auto f = [&](double z) {
                return func(x, y, z);
            };
            return qgaus(f, z1(x, y), z2(x, y));
        };
        return qgaus(g, y1(x), y2(x));
    };
    return qgaus(h, x1, x2);
}
