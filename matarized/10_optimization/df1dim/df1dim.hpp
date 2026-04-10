#pragma once
#include <matar.h>

using namespace mtr;

// 1-D directional-derivative evaluator for gradient-based line searches.
// Replaces the Fortran COMMON-block pattern. Stores base point p, direction xi,
// n-dimensional gradient function dfunc. Evaluating operator()(x) returns the
// directional derivative dot(grad f(p + x*xi), xi).
//
// Preferred usage: construct a lambda inside dlinmin. This struct is provided
// for standalone use or when a named type is needed.
template <typename DFunc>
struct DF1Dim {
    const double* p;
    const double* xi;
    int           n;
    DFunc         dfunc;

    DF1Dim(const double* p_, const double* xi_, int n_, DFunc dfunc_)
        : p(p_), xi(xi_), n(n_), dfunc(dfunc_) {}

    double operator()(double x) const
    {
        constexpr int NMAX = 50;
        double xt[NMAX], df[NMAX];
        for (int j = 0; j < n; j++)
            xt[j] = p[j] + x * xi[j];
        dfunc(xt, df);
        double val = 0.0;
        for (int j = 0; j < n; j++)
            val += df[j] * xi[j];
        return val;
    }
};

// Convenience factory
template <typename DFunc>
inline DF1Dim<DFunc> make_df1dim(const double* p, const double* xi,
                                  int n, DFunc dfunc)
{
    return DF1Dim<DFunc>(p, xi, n, dfunc);
}
