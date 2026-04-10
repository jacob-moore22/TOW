// Singular Value Decomposition (Numerical Recipes SVDCMP).
//
// Given a(1..m, 1..n) with m >= n, computes its SVD:  A = U * W * V^T
// On return:
//   a  is replaced by U  (m x n)
//   w  contains the n singular values (1..n)
//   v  contains V (not V^T) as an n x n matrix (1..n, 1..n)

#pragma once
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <matar.h>

using namespace mtr;

namespace detail {
inline double pythag(double a, double b)
{
    double absa = std::fabs(a);
    double absb = std::fabs(b);
    if (absa > absb) {
        double r = absb / absa;
        return absa * std::sqrt(1.0 + r * r);
    }
    if (absb == 0.0) return 0.0;
    double r = absa / absb;
    return absb * std::sqrt(1.0 + r * r);
}
} // namespace detail

inline void svdcmp(DFMatrixKokkos<double>& a, int m, int n,
                   DFMatrixKokkos<double>& w, DFMatrixKokkos<double>& v)
{
    double rv1[n];
    double g = 0.0, scale = 0.0, anorm = 0.0;

    // Householder reduction to bidiagonal form
    for (int i = 1; i <= n; i++) {
        int l = i + 1;
        rv1[i - 1] = scale * g;
        g = 0.0;
        double s = 0.0;
        scale = 0.0;

        if (i <= m) {
            for (int k = i; k <= m; k++)
                scale += std::fabs(a.host(k, i));
            if (scale != 0.0) {
                for (int k = i; k <= m; k++) {
                    a.host(k, i) /= scale;
                    s += a.host(k, i) * a.host(k, i);
                }
                double f = a.host(i, i);
                g = -std::copysign(std::sqrt(s), f);
                double h = f * g - s;
                a.host(i, i) = f - g;
                if (i != n) {
                    for (int j = l; j <= n; j++) {
                        s = 0.0;
                        for (int k = i; k <= m; k++)
                            s += a.host(k, i) * a.host(k, j);
                        f = s / h;
                        for (int k = i; k <= m; k++)
                            a.host(k, j) += f * a.host(k, i);
                    }
                }
                for (int k = i; k <= m; k++)
                    a.host(k, i) *= scale;
            }
        }

        w.host(i) = scale * g;
        g = 0.0;
        s = 0.0;
        scale = 0.0;

        if (i <= m && i != n) {
            for (int k = l; k <= n; k++)
                scale += std::fabs(a.host(i, k));
            if (scale != 0.0) {
                for (int k = l; k <= n; k++) {
                    a.host(i, k) /= scale;
                    s += a.host(i, k) * a.host(i, k);
                }
                double f = a.host(i, l);
                g = -std::copysign(std::sqrt(s), f);
                double h = f * g - s;
                a.host(i, l) = f - g;
                for (int k = l; k <= n; k++)
                    rv1[k - 1] = a.host(i, k) / h;
                if (i != m) {
                    for (int j = l; j <= m; j++) {
                        s = 0.0;
                        for (int k = l; k <= n; k++)
                            s += a.host(j, k) * a.host(i, k);
                        for (int k = l; k <= n; k++)
                            a.host(j, k) += s * rv1[k - 1];
                    }
                }
                for (int k = l; k <= n; k++)
                    a.host(i, k) *= scale;
            }
        }

        anorm = std::fmax(anorm, std::fabs(w.host(i)) + std::fabs(rv1[i - 1]));
    }

    // Accumulation of right-hand transformations
    for (int i = n; i >= 1; i--) {
        int l = i + 1;
        if (i < n) {
            if (g != 0.0) {
                for (int j = l; j <= n; j++)
                    v.host(j, i) = (a.host(i, j) / a.host(i, l)) / g;
                for (int j = l; j <= n; j++) {
                    double s = 0.0;
                    for (int k = l; k <= n; k++)
                        s += a.host(i, k) * v.host(k, j);
                    for (int k = l; k <= n; k++)
                        v.host(k, j) += s * v.host(k, i);
                }
            }
            for (int j = l; j <= n; j++) {
                v.host(i, j) = 0.0;
                v.host(j, i) = 0.0;
            }
        }
        v.host(i, i) = 1.0;
        g = rv1[i - 1];
    }

    // Accumulation of left-hand transformations
    for (int i = n; i >= 1; i--) {
        int l = i + 1;
        g = w.host(i);
        if (i < n) {
            for (int j = l; j <= n; j++)
                a.host(i, j) = 0.0;
        }
        if (g != 0.0) {
            g = 1.0 / g;
            if (i != n) {
                for (int j = l; j <= n; j++) {
                    double s = 0.0;
                    for (int k = l; k <= m; k++)
                        s += a.host(k, i) * a.host(k, j);
                    double f = (s / a.host(i, i)) * g;
                    for (int k = i; k <= m; k++)
                        a.host(k, j) += f * a.host(k, i);
                }
            }
            for (int j = i; j <= m; j++)
                a.host(j, i) *= g;
        } else {
            for (int j = i; j <= m; j++)
                a.host(j, i) = 0.0;
        }
        a.host(i, i) += 1.0;
    }

    // Diagonalization of the bidiagonal form
    for (int k = n; k >= 1; k--) {
        for (int its = 1; its <= 30; its++) {
            bool flag = true;
            int l, nm = 0;
            for (l = k; l >= 1; l--) {
                nm = l - 1;
                if (std::fabs(rv1[l - 1]) + anorm == anorm) { flag = false; break; }
                if (nm >= 1 && std::fabs(w.host(nm)) + anorm == anorm) break;
            }

            if (flag) {
                double c = 0.0, s = 1.0;
                for (int i = l; i <= k; i++) {
                    double f = s * rv1[i - 1];
                    rv1[i - 1] *= c;
                    if (std::fabs(f) + anorm == anorm) break;
                    g = w.host(i);
                    double h = detail::pythag(f, g);
                    w.host(i) = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (int j = 1; j <= m; j++) {
                        double y = a.host(j, nm);
                        double z = a.host(j, i);
                        a.host(j, nm) = y * c + z * s;
                        a.host(j, i)  = -y * s + z * c;
                    }
                }
            }

            double z = w.host(k);
            if (l == k) {
                if (z < 0.0) {
                    w.host(k) = -z;
                    for (int j = 1; j <= n; j++)
                        v.host(j, k) = -v.host(j, k);
                }
                break;
            }

            if (its == 30) {
                std::fprintf(stderr, "svdcmp: no convergence in 30 iterations\n");
                return;
            }

            double x = w.host(l);
            nm = k - 1;
            double y = w.host(nm);
            g = rv1[nm - 1];
            double h = rv1[k - 1];
            double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = detail::pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + std::copysign(g, f))) - h)) / x;

            double c = 1.0, s = 1.0;
            for (int j = l; j <= nm; j++) {
                int i = j + 1;
                g = rv1[i - 1];
                y = w.host(i);
                h = s * g;
                g = c * g;
                z = detail::pythag(f, h);
                rv1[j - 1] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = -x * s + g * c;
                h = y * s;
                y = y * c;
                for (int jj = 1; jj <= n; jj++) {
                    x = v.host(jj, j);
                    z = v.host(jj, i);
                    v.host(jj, j) = x * c + z * s;
                    v.host(jj, i) = -x * s + z * c;
                }
                z = detail::pythag(f, h);
                w.host(j) = z;
                if (z != 0.0) {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = -s * g + c * y;
                for (int jj = 1; jj <= m; jj++) {
                    y = a.host(jj, j);
                    z = a.host(jj, i);
                    a.host(jj, j) = y * c + z * s;
                    a.host(jj, i) = -y * s + z * c;
                }
            }
            rv1[l - 1] = 0.0;
            rv1[k - 1] = f;
            w.host(k) = x;
        }
    }
}
