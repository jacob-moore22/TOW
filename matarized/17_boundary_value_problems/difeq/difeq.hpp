#pragma once
#include <matar.h>

using namespace mtr;

// Generic finite-difference equation interface for solvde
// (Numerical Recipes DIFEQ).
//
// This is problem-specific: each BVP supplies its own Difeq callable.
// The callable must have signature:
//   void difeq(int k, int k1, int k2, int jsf,
//              int is1, int isf, int* indexv, int ne,
//              DFMatrixKokkos<double>& s, int nsi, int nsj,
//              DFMatrixKokkos<double>& y, int nyj, int nyk)
//
// For the spheroidal harmonics test problem, see sfroid.hpp.

// Spheroidal-wave difeq functor used by sfroid.
// Captures problem parameters (x grid, h, mm, n, c2, anorm) by reference.
struct SpheroidalDifeq {
    double* x;
    double  h;
    int     mm;
    int     n;
    double  c2;
    double  anorm;
    int     m_grid;

    void operator()(int k, int k1, int k2, int jsf,
                    int is1, int isf, int* indexv, int ne,
                    DFMatrixKokkos<double>& s, int nsi, int nsj,
                    DFMatrixKokkos<double>& y, int nyj, int nyk) const
    {
        for (int i = 1; i <= nsi; i++)
            for (int j = 1; j <= nsj; j++)
                s.host(i, j) = 0.0;

        if (k == k1) {
            if ((n + mm) % 2 == 1) {
                s.host(3, 3 + indexv[1]) = 1.0;
                s.host(3, jsf) = y.host(1, 1);
            } else {
                s.host(3, 3 + indexv[2]) = 1.0;
                s.host(3, jsf) = y.host(2, 1);
            }
        } else if (k > k2) {
            s.host(1, 3 + indexv[1]) = -(y.host(3, m_grid) - c2) / (2.0 * (mm + 1.0));
            s.host(1, 3 + indexv[2]) = 1.0;
            s.host(1, 3 + indexv[3]) = -y.host(1, m_grid) / (2.0 * (mm + 1.0));
            s.host(1, jsf) = y.host(2, m_grid)
                           - (y.host(3, m_grid) - c2) * y.host(1, m_grid)
                             / (2.0 * (mm + 1.0));
            s.host(2, 3 + indexv[1]) = 1.0;
            s.host(2, jsf) = y.host(1, m_grid) - anorm;
        } else {
            s.host(1, indexv[1])     = -1.0;
            s.host(1, indexv[2])     = -0.5 * h;
            s.host(1, 3 + indexv[1]) =  1.0;
            s.host(1, 3 + indexv[2]) = -0.5 * h;

            double xmid = 0.5 * (x[k] + x[k - 1]);
            double temp  = h / (1.0 - xmid * xmid);
            double temp2 = 0.5 * (y.host(3, k) + y.host(3, k - 1))
                         - c2 * xmid * xmid;

            s.host(2, indexv[1])     =  0.5 * temp * temp2;
            s.host(2, indexv[2])     = -1.0 - 0.5 * temp * (mm + 1.0) * (x[k] + x[k - 1]);
            s.host(2, indexv[3])     =  0.25 * temp * (y.host(1, k) + y.host(1, k - 1));
            s.host(2, 3 + indexv[1]) =  s.host(2, indexv[1]);
            s.host(2, 3 + indexv[2]) =  2.0 + s.host(2, indexv[2]);
            s.host(2, 3 + indexv[3]) =  s.host(2, indexv[3]);

            s.host(3, indexv[3])     = -1.0;
            s.host(3, 3 + indexv[3]) =  1.0;

            s.host(1, jsf) = y.host(1, k) - y.host(1, k - 1)
                           - 0.5 * h * (y.host(2, k) + y.host(2, k - 1));
            s.host(2, jsf) = y.host(2, k) - y.host(2, k - 1)
                           - temp * ((x[k] + x[k - 1]) * 0.5 * (mm + 1.0)
                                     * (y.host(2, k) + y.host(2, k - 1))
                                     - temp2 * 0.5
                                     * (y.host(1, k) + y.host(1, k - 1)));
            s.host(3, jsf) = y.host(3, k) - y.host(3, k - 1);
        }
    }
};
