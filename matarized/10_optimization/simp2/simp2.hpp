#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Helper for the simplex LP solver. Determines the pivot row ip for column kp
// using the minimum-ratio test with Bland-style degeneracy breaking.
//
// a:    tableau, row-major, dimensions (mp rows) x (np cols).
// m, n: number of constraints / variables.
// np:   column stride of the tableau.
// l2:   array of row indices to consider, length nl2.
// nl2:  number of entries in l2.
// ip:   on return, pivot row (0 if unbounded).
// kp:   pivot column (0-based).
// q1:   on return, the minimum ratio.
inline void simp2(const double* a, int m, int n, int np,
                  const int* l2, int nl2, int& ip, int kp, double& q1)
{
    constexpr double EPS = 1.0e-6;
    ip = 0;
    if (nl2 < 1) return;

    int first = -1;
    for (int i = 0; i < nl2; i++) {
        if (a[l2[i] * np + kp] < -EPS) {
            first = i;
            break;
        }
    }
    if (first < 0) return;

    q1 = -a[l2[first] * np + 0] / a[l2[first] * np + kp];
    ip = l2[first];

    for (int i = first + 1; i < nl2; i++) {
        int ii = l2[i];
        if (a[ii * np + kp] < -EPS) {
            double q = -a[ii * np + 0] / a[ii * np + kp];
            if (q < q1) {
                ip = ii;
                q1 = q;
            } else if (q == q1) {
                // Bland's rule tie-breaking
                double qp = 0.0, q0 = 0.0;
                for (int k = 1; k <= n; k++) {
                    qp = -a[ip * np + k] / a[ip * np + kp];
                    q0 = -a[ii * np + k] / a[ii * np + kp];
                    if (q0 != qp) break;
                }
                if (q0 < qp) ip = ii;
            }
        }
    }
}
