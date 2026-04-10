#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Shooting to a fitting point (Numerical Recipes SHOOTF).
//
// Template params replace Fortran EXTERNAL routines:
//   Load1:  void load1(double x1, DFMatrixKokkos<double>& v1, DFMatrixKokkos<double>& y)
//   Load2:  void load2(double x2, DFMatrixKokkos<double>& v2, DFMatrixKokkos<double>& y)
//   Score:  void score(double xf, DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& f)
//   Derivs: void derivs(double x, DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx)
//   Stepper: rkqc-compatible
//
// v1(1..n2): free params at x1;  v2(1..n1): free params at x2
// Requires: odeint.hpp, ludcmp.hpp, lubksb.hpp
template <typename Load1, typename Load2, typename Score,
          typename Derivs, typename Stepper>
inline void shootf(int nvar,
                   DFMatrixKokkos<double>& v1, DFMatrixKokkos<double>& v2,
                   DFMatrixKokkos<double>& delv1, DFMatrixKokkos<double>& delv2,
                   int n1, int n2,
                   double x1, double x2, double xf,
                   double eps, double h1, double hmin,
                   DFMatrixKokkos<double>& f,
                   DFMatrixKokkos<double>& dv1, DFMatrixKokkos<double>& dv2,
                   Load1 load1, Load2 load2, Score score,
                   Derivs derivs, Stepper stepper)
{
    DFMatrixKokkos<double> y(nvar);
    DFMatrixKokkos<double> f1(nvar);
    DFMatrixKokkos<double> f2(nvar);
    DFMatrixKokkos<double> dfdv(nvar, nvar);

    OdeintPath path;
    int nok, nbad;

    // Integrate from x1 to xf
    load1(x1, v1, y);
    odeint(y, nvar, x1, xf, eps, h1, hmin, nok, nbad, derivs, stepper, path);
    score(xf, y, f1);

    // Integrate from x2 to xf
    load2(x2, v2, y);
    odeint(y, nvar, x2, xf, eps, h1, hmin, nok, nbad, derivs, stepper, path);
    score(xf, y, f2);

    // Jacobian from v1 perturbations
    int j = 0;
    for (int iv = 1; iv <= n2; iv++) {
        j++;
        double sav = v1.host(iv);
        v1.host(iv) += delv1.host(iv);
        load1(x1, v1, y);
        odeint(y, nvar, x1, xf, eps, h1, hmin, nok, nbad, derivs, stepper, path);
        score(xf, y, f);
        for (int i = 1; i <= nvar; i++)
            dfdv.host(i, j) = (f.host(i) - f1.host(i)) / delv1.host(iv);
        v1.host(iv) = sav;
    }

    // Jacobian from v2 perturbations (sign reversed)
    for (int iv = 1; iv <= n1; iv++) {
        j++;
        double sav = v2.host(iv);
        v2.host(iv) += delv2.host(iv);
        load2(x2, v2, y);
        odeint(y, nvar, x2, xf, eps, h1, hmin, nok, nbad, derivs, stepper, path);
        score(xf, y, f);
        for (int i = 1; i <= nvar; i++)
            dfdv.host(i, j) = (f2.host(i) - f.host(i)) / delv2.host(iv);
        v2.host(iv) = sav;
    }

    // RHS = -(f1 - f2)
    for (int i = 1; i <= nvar; i++) {
        f.host(i)  = f1.host(i) - f2.host(i);
        f1.host(i) = -f.host(i);
    }

    // Solve for corrections via LU
    DFMatrixKokkos<double> a_mat(nvar, nvar);
    DFMatrixKokkos<int>    indx_mat(nvar);
    DFMatrixKokkos<double> b_vec(nvar);

    for (int i = 1; i <= nvar; i++) {
        for (int jj = 1; jj <= nvar; jj++)
            a_mat.host(i, jj) = dfdv.host(i, jj);
        b_vec.host(i) = f1.host(i);
    }

    double det;
    ludcmp(a_mat, nvar, indx_mat, det);
    lubksb(a_mat, nvar, indx_mat, b_vec);

    // Distribute corrections
    j = 0;
    for (int iv = 1; iv <= n2; iv++) {
        j++;
        v1.host(iv)  += b_vec.host(j);
        dv1.host(iv)  = b_vec.host(j);
    }
    for (int iv = 1; iv <= n1; iv++) {
        j++;
        v2.host(iv)  += b_vec.host(j);
        dv2.host(iv)  = b_vec.host(j);
    }
}
