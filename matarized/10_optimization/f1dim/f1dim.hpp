#pragma once
#include <matar.h>

using namespace mtr;

// 1-D function evaluator for line searches. Replaces the Fortran COMMON-block
// pattern: stores base point p, direction xi, and the n-dimensional objective
// func. Evaluating operator()(x) returns func(p + x * xi).
//
// Preferred usage: construct a lambda inside linmin. This struct is provided
// for standalone use or when a named type is needed.
template <typename Func>
struct F1Dim {
    const double* p;
    const double* xi;
    int           n;
    Func          func;

    F1Dim(const double* p_, const double* xi_, int n_, Func func_)
        : p(p_), xi(xi_), n(n_), func(func_) {}

    double operator()(double x) const
    {
        constexpr int NMAX = 50;
        double xt[NMAX];
        for (int j = 0; j < n; j++)
            xt[j] = p[j] + x * xi[j];
        return func(xt);
    }
};

// Convenience factory
template <typename Func>
inline F1Dim<Func> make_f1dim(const double* p, const double* xi, int n, Func func)
{
    return F1Dim<Func>(p, xi, n, func);
}
