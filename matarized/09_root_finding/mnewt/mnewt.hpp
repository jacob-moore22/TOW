#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// LU decomposition (Crout's algorithm with partial pivoting).
// a is n x n (0-indexed), indx[0..n-1] receives pivot indices, d is +/-1.
inline void mnewt_ludcmp(double** a, int n, int* indx, double& d)
{
    constexpr double TINY = 1.0e-20;
    double vv[64];

    d = 1.0;
    for (int i = 0; i < n; i++) {
        double big = 0.0;
        for (int j = 0; j < n; j++) {
            double temp = std::fabs(a[i][j]);
            if (temp > big) big = temp;
        }
        if (big == 0.0) {
            std::fprintf(stderr, "mnewt_ludcmp: singular matrix\n");
            return;
        }
        vv[i] = 1.0 / big;
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            double sum = a[i][j];
            for (int k = 0; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        double big = 0.0;
        int imax = j;
        for (int i = j; i < n; i++) {
            double sum = a[i][j];
            for (int k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            double dum = vv[i] * std::fabs(sum);
            if (dum >= big) {
                big  = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (int k = 0; k < n; k++) {
                double dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k]    = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n - 1) {
            double dum = 1.0 / a[j][j];
            for (int i = j + 1; i < n; i++)
                a[i][j] *= dum;
        }
    }
}

// LU back-substitution. a is the LU-decomposed matrix, indx from ludcmp,
// b[0..n-1] is the RHS on entry, solution on exit.
inline void mnewt_lubksb(double** a, int n, const int* indx, double* b)
{
    int ii = -1;
    for (int i = 0; i < n; i++) {
        int ll = indx[i];
        double sum = b[ll];
        b[ll] = b[i];
        if (ii >= 0) {
            for (int j = ii; j < i; j++)
                sum -= a[i][j] * b[j];
        } else if (sum != 0.0) {
            ii = i;
        }
        b[i] = sum;
    }
    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

// Multi-dimensional Newton's method. usrfun(x, alpha, beta) computes the
// Jacobian alpha[i][j] = d(f_i)/d(x_j) and the function values beta[i] = -f_i.
// x[0..n-1] is the initial guess, refined in place over ntrial iterations
// until ||f|| < tolf or ||dx|| < tolx.
template<typename UserFun>
inline void mnewt(int ntrial, double* x, int n, double tolx, double tolf,
                  UserFun usrfun)
{
    constexpr int NP = 64;
    double  alpha_data[NP][NP];
    double* alpha[NP];
    for (int i = 0; i < n; i++) alpha[i] = alpha_data[i];
    double beta[NP];
    int    indx[NP];

    for (int k = 0; k < ntrial; k++) {
        usrfun(x, alpha, beta);

        double errf = 0.0;
        for (int i = 0; i < n; i++)
            errf += std::fabs(beta[i]);
        if (errf <= tolf) return;

        double d;
        mnewt_ludcmp(alpha, n, indx, d);
        mnewt_lubksb(alpha, n, indx, beta);

        double errx = 0.0;
        for (int i = 0; i < n; i++) {
            errx += std::fabs(beta[i]);
            x[i] += beta[i];
        }
        if (errx <= tolx) return;
    }
}
