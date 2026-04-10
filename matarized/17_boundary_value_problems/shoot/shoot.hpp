#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Shooting method for two-point BVP (Numerical Recipes SHOOT).
//
// Template params replace Fortran EXTERNAL routines:
//   Load:   void load(double x1, DFMatrixKokkos<double>& v, DFMatrixKokkos<double>& y)
//   Score:  void score(double x2, DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& f)
//   Derivs: void derivs(double x, DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx)
//   Stepper: rkqc-compatible
//
// v(1..n2): adjustable parameters at x1
// delv(1..n2): step sizes for Jacobian
// f(1..n2), dv(1..n2): workspace/output corrections
//
// Requires: odeint.hpp (with OdeintPath), ludcmp.hpp, lubksb.hpp
template <typename Load, typename Score, typename Derivs, typename Stepper>
inline void shoot(int nvar, DFMatrixKokkos<double>& v,
                  DFMatrixKokkos<double>& delv, int n2,
                  double x1, double x2, double eps, double h1, double hmin,
                  DFMatrixKokkos<double>& f, DFMatrixKokkos<double>& dv,
                  Load load_fn, Score score_fn, Derivs derivs, Stepper stepper)
{
    DFMatrixKokkos<double> y(nvar);
    DFMatrixKokkos<double> dfdv(n2, n2);

    OdeintPath path;
    int nok, nbad;

    load_fn(x1, v, y);
    odeint(y, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, stepper, path);
    score_fn(x2, y, f);

    for (int iv = 1; iv <= n2; iv++) {
        double sav = v.host(iv);
        v.host(iv) += delv.host(iv);
        load_fn(x1, v, y);
        odeint(y, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, stepper, path);
        score_fn(x2, y, dv);
        for (int i = 1; i <= n2; i++)
            dfdv.host(i, iv) = (dv.host(i) - f.host(i)) / delv.host(iv);
        v.host(iv) = sav;
    }

    for (int iv = 1; iv <= n2; iv++)
        dv.host(iv) = -f.host(iv);

    DFMatrixKokkos<double> a_mat(n2, n2);
    DFMatrixKokkos<int>    indx_mat(n2);
    DFMatrixKokkos<double> b_vec(n2);

    for (int i = 1; i <= n2; i++) {
        for (int j = 1; j <= n2; j++)
            a_mat.host(i, j) = dfdv.host(i, j);
        b_vec.host(i) = dv.host(i);
    }

    double det;
    ludcmp(a_mat, n2, indx_mat, det);
    lubksb(a_mat, n2, indx_mat, b_vec);

    for (int iv = 1; iv <= n2; iv++)
        v.host(iv) += b_vec.host(iv);
}
