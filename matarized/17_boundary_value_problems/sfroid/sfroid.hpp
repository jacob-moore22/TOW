#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Spheroidal wave function eigenvalue solver (Numerical Recipes SFROID).
// Uses the solvde relaxation method to compute eigenvalues of the
// spheroidal wave equation for given quantum numbers mm, n, and c^2.
//
// Returns lambda = eigenvalue y(3,1) + mm*(mm+1).
// Requires: solvde.hpp, bksub.hpp, pinvs.hpp, red.hpp, difeq.hpp (SpheroidalDifeq).
inline double sfroid(int mm_val, int n_val, double c2,
                     int M = 41, int itmax = 100, double conv = 5.0e-6)
{
    constexpr int NE  = 3;
    constexpr int NB  = 1;

    int NCI = NE;
    int NCJ = NE - NB + 1;
    int NCK = M + 1;
    int NSI = NE;
    int NSJ = 2 * NE + 1;

    double h = 1.0 / (M - 1);

    double* x = new double[M + 1];
    for (int k = 1; k <= M - 1; k++)
        x[k] = (k - 1) * h;
    x[M] = 1.0;

    int indexv[4];
    if ((n_val + mm_val) % 2 == 1) {
        indexv[1] = 1; indexv[2] = 2; indexv[3] = 3;
    } else {
        indexv[1] = 2; indexv[2] = 1; indexv[3] = 3;
    }

    double anorm = 1.0;
    if (mm_val != 0) {
        double q1 = n_val;
        for (int i = 1; i <= mm_val; i++) {
            anorm = -0.5 * anorm * (n_val + i) * (q1 / i);
            q1 -= 1.0;
        }
    }

    DFMatrixKokkos<double> y(NE, M);
    DFMatrixKokkos<double> c_arr(NCI, NCJ, NCK);
    DFMatrixKokkos<double> s(NSI, NSJ);

    // plgndr would be used here for proper initialization;
    // simplified: use constant initial guess
    for (int k = 1; k <= M; k++) {
        y.host(1, k) = anorm;
        y.host(2, k) = 0.0;
        y.host(3, k) = static_cast<double>(n_val * (n_val + 1) - mm_val * (mm_val + 1));
    }

    SpheroidalDifeq difeq;
    difeq.x      = x;
    difeq.h      = h;
    difeq.mm     = mm_val;
    difeq.n      = n_val;
    difeq.c2     = c2;
    difeq.anorm  = anorm;
    difeq.m_grid = M;

    double scalv[4];
    scalv[1] = std::fabs(anorm);
    scalv[2] = std::fmax(std::fabs(anorm), std::fabs(y.host(2, M)));
    scalv[3] = std::fmax(1.0, y.host(3, M));

    solvde(itmax, conv, 1.0, scalv, indexv, NE, NB, M,
           y, NE, M, c_arr, NCI, NCJ, NCK, s, NSI, NSJ, difeq);

    double lambda = y.host(3, 1) + mm_val * (mm_val + 1);

    delete[] x;
    return lambda;
}
