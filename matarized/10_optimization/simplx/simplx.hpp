#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "../simp1/simp1.hpp"
#include "../simp2/simp2.hpp"
#include "../simp3/simp3.hpp"

using namespace mtr;

// Simplex method for linear programming.
//
// Maximize z = sum_j c_j x_j subject to:
//   sum_j a_ij x_j <= b_i   for i = 1..m1
//   sum_j a_ij x_j >= b_i   for i = m1+1..m1+m2
//   sum_j a_ij x_j  = b_i   for i = m1+m2+1..m
//
// All b_i must be >= 0.  m = m1 + m2 + m3.
//
// Tableau layout (row-major, 0-based C++ storage, but using NR 1-based
// indexing conventions internally):
//   a[(m+2) rows x (n+1) cols]
//   Row 0:     [z0, c_1, ..., c_n]           objective
//   Rows 1..m: [b_i, -a_i1, ..., -a_in]      constraints (note negated a_ij)
//   Row m+1:   [0, 0, ..., 0]                auxiliary (workspace)
//
// np = n+1 is the column stride.
//
// On return:
//   icase =  0 : optimal solution found
//   icase =  1 : unbounded
//   icase = -1 : no feasible solution
//
//   izrov[0..n-1]:   1-based indices of non-basic variables
//   iposv[0..m-1]:   1-based indices of basic variables
//
//   Optimal z is in a[0].  For each i, if iposv[i] is a 1-based variable
//   index k (1 <= k <= n), then x_k = a[(i+1)*np + 0].
inline void simplx(double* a, int m, int n, int np,
                   int m1, int m2, int m3,
                   int& icase, int* izrov, int* iposv)
{
    constexpr int    MMAX = 100;
    constexpr double EPS  = 1.0e-6;

    int l1[MMAX], l2[MMAX], l3[MMAX];

    if (m != m1 + m2 + m3) {
        std::fprintf(stderr, "simplx: bad constraint counts\n");
        return;
    }

    int nl1 = n;
    for (int k = 0; k < n; k++) {
        l1[k]    = k + 1;
        izrov[k] = k + 1;
    }
    int nl2 = m;
    for (int i = 0; i < m; i++) {
        if (a[(i + 1) * np + 0] < 0.0) {
            std::fprintf(stderr, "simplx: bad input tableau (negative RHS)\n");
            return;
        }
        l2[i]    = i + 1;
        iposv[i] = n + i + 1;
    }
    for (int i = 0; i < m2; i++) l3[i] = 1;

    int ir = 0;
    int kp = 0, ip = 0;
    double bmax = 0.0, q1 = 0.0;

    if (m2 + m3 == 0) goto label_30;

    // Phase 1: auxiliary objective
    ir = 1;
    for (int k = 0; k <= n; k++) {
        q1 = 0.0;
        for (int i = m1 + 1; i <= m; i++)
            q1 += a[i * np + k];
        a[(m + 1) * np + k] = -q1;
    }

label_10:
    simp1(a, np, m + 1, l1, nl1, 0, kp, bmax);
    if (bmax <= EPS && a[(m + 1) * np + 0] < -EPS) {
        icase = -1;
        return;
    } else if (bmax <= EPS && a[(m + 1) * np + 0] <= EPS) {
        // Feasible solution found for auxiliary problem
        {
            int m12 = m1 + m2 + 1;
            if (m12 <= m) {
                for (int iip = m12; iip <= m; iip++) {
                    if (iposv[iip - 1] == iip + n) {
                        simp1(a, np, iip, l1, nl1, 1, kp, bmax);
                        if (bmax > 0.0) { ip = iip; goto label_1; }
                    }
                }
            }
        }
        ir = 0;
        {
            int m12 = m1 + m2;
            if (m1 + 1 <= m12) {
                for (int i = m1 + 1; i <= m12; i++) {
                    if (l3[i - m1 - 1] == 1) {
                        for (int k = 0; k <= n; k++)
                            a[i * np + k] = -a[i * np + k];
                    }
                }
            }
        }
        goto label_30;
    }

    simp2(a, m, n, np, l2, nl2, ip, kp, q1);
    if (ip == 0) {
        icase = -1;
        return;
    }

label_1:
    simp3(a, np, m + 1, n, ip, kp);

    if (iposv[ip - 1] >= n + m1 + m2 + 1) {
        // Auxiliary variable leaving: remove kp from l1
        int ki = -1;
        for (int k = 0; k < nl1; k++) {
            if (l1[k] == kp) { ki = k; break; }
        }
        nl1--;
        for (int is = ki; is < nl1; is++)
            l1[is] = l1[is + 1];
    } else {
        if (iposv[ip - 1] < n + m1 + 1) goto label_20;
        int kh = iposv[ip - 1] - m1 - n;
        if (l3[kh - 1] == 0) goto label_20;
        l3[kh - 1] = 0;
    }

    a[(m + 1) * np + kp] += 1.0;
    for (int i = 0; i <= m + 1; i++)
        a[i * np + kp] = -a[i * np + kp];

label_20:
    {
        int is        = izrov[kp - 1];
        izrov[kp - 1] = iposv[ip - 1];
        iposv[ip - 1] = is;
    }
    if (ir != 0) goto label_10;

label_30:
    simp1(a, np, 0, l1, nl1, 0, kp, bmax);
    if (bmax <= 0.0) {
        icase = 0;
        return;
    }
    simp2(a, m, n, np, l2, nl2, ip, kp, q1);
    if (ip == 0) {
        icase = 1;
        return;
    }
    simp3(a, np, m, n, ip, kp);
    goto label_20;
}
