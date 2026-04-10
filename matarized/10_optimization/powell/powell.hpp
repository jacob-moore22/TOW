#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "../linmin/linmin.hpp"

using namespace mtr;

// Powell's method of conjugate directions (derivative-free).
//
// p:    starting point (length n), replaced by the minimum on return.
// xi:   initial direction set, n×n stored row-major (xi[i*n+j] is direction i,
//       component j). Updated on return.
// n:    number of dimensions.
// ftol: fractional convergence tolerance.
// iter: on return, number of iterations performed.
// fret: on return, minimum function value.
// func: objective f(const double* x) -> double.
template <typename Func>
inline void powell(double* p, double* xi, int n, double ftol,
                   int& iter, double& fret, Func func)
{
    constexpr int NMAX  = 50;
    constexpr int ITMAX = 200;

    double pt[NMAX], ptt[NMAX], xit[NMAX];

    fret = func(p);
    for (int j = 0; j < n; j++) pt[j] = p[j];

    iter = 0;
    for (;;) {
        iter++;
        double fp = fret;
        int ibig = 0;
        double del = 0.0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) xit[j] = xi[i * n + j];
            double fptt = fret;
            linmin(p, xit, n, fret, func);
            if (std::fabs(fptt - fret) > del) {
                del  = std::fabs(fptt - fret);
                ibig = i;
            }
        }

        if (2.0 * std::fabs(fp - fret) <= ftol * (std::fabs(fp) + std::fabs(fret)))
            return;

        if (iter >= ITMAX) {
            std::fprintf(stderr, "powell: exceeded maximum iterations\n");
            return;
        }

        for (int j = 0; j < n; j++) {
            ptt[j] = 2.0 * p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j]  = p[j];
        }

        double fptt = func(ptt);
        if (fptt >= fp) continue;

        double t = 2.0 * (fp - 2.0 * fret + fptt) *
                   (fp - fret - del) * (fp - fret - del) -
                   del * (fp - fptt) * (fp - fptt);
        if (t >= 0.0) continue;

        linmin(p, xit, n, fret, func);
        for (int j = 0; j < n; j++)
            xi[ibig * n + j] = xit[j];
    }
}
